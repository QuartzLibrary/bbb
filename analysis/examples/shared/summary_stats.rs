use std::{cmp::Ordering, str::FromStr};

use analysis::util::InspectEvery;
use futures::{StreamExt, stream};
use hail::contig::GRCh38Contig;
use pan_ukbb::{PhenotypeManifestEntry, SummaryStats};
use tokio::task::spawn_blocking;
use utile::{
    cache::FsCacheEntry,
    resource::{RawResource, RawResourceExt},
};

type FastaReader = biocore::fasta::IndexedFastaReader<std::io::BufReader<std::fs::File>>;

pub async fn merged_summary_stats(
    phenotypes: Vec<PhenotypeManifestEntry>,
    lock: &'static tokio::sync::Mutex<()>,
    liftover: &'static liftover::LiftoverIndexed,
) -> impl Iterator<Item = (String, SummaryStats<GRCh38Contig>)> {
    let all_summary_stats: Vec<_> = stream::iter(phenotypes.clone())
        .map(|p| async move {
            let filename = p.filename.clone();
            (filename, summary_stats(p.clone(), liftover, lock).await)
        })
        .buffered(100)
        .collect()
        .await;

    let all_summary_stats = all_summary_stats
        .into_iter()
        .map(|(filename, stats)| stats.map(move |s| (filename.clone(), s)));

    itertools::kmerge_by(
        all_summary_stats,
        |a: &(String, SummaryStats<GRCh38Contig>), b: &(String, SummaryStats<GRCh38Contig>)| {
            let (_, a) = a;
            let (_, b) = b;
            cmp_summary_stats(a, b) == Ordering::Less
        },
    )
}

pub async fn summary_stats(
    phenotype: PhenotypeManifestEntry,
    liftover: &'static liftover::LiftoverIndexed,
    lock: &'static tokio::sync::Mutex<()>,
) -> impl Iterator<Item = SummaryStats<GRCh38Contig>> {
    let summary_stats_entry = super::summary_stats_temp_path(&phenotype);

    if !summary_stats_entry.try_exists().unwrap() {
        load_and_cache_summary_stats(phenotype.clone(), &summary_stats_entry, liftover, lock).await;
    }

    summary_stats_entry
        .clone()
        .decompressed_with(utile::resource::Compression::Gzip)
        .read_json_lines()
        .unwrap()
        .map(|s| s.unwrap())
}

async fn load_and_cache_summary_stats(
    phenotype: PhenotypeManifestEntry,
    entry: &FsCacheEntry,
    liftover: &'static liftover::LiftoverIndexed,
    lock: &'static tokio::sync::Mutex<()>,
) {
    log::info!("[SummaryStats][{}] Starting", phenotype.filename);

    let mut grch38: FastaReader = hail::load_grch38_reference_genome().await.unwrap();
    let mut grch37 = hail::load_grch37_reference_genome().await.unwrap();
    log::info!("[SummaryStats][{}] Loaded refs", phenotype.filename);

    let summary_stats = SummaryStats::<hail::contig::GRCh37Contig>::load(
        phenotype
            .summary_stats_resource()
            .log_progress()
            .with_fs_cache(&super::DATA_FOLDER)
            .ensure_cached_async()
            .await
            .unwrap()
            .decompressed()
            .buffered(),
    )
    .unwrap()
    .map(|s| s.unwrap());

    log::info!("[SummaryStats][{}] Cached", phenotype.filename);

    let lock = lock.lock().await;

    spawn_blocking({
        let entry = entry.clone();
        move || {
            std::thread::spawn(move || {
                let summary_stats =
                    load_summary_stats(summary_stats, &mut grch38, &mut grch37, liftover);

                entry
                    .write_file(
                        utile::resource::iter::IterToJsonLinesResource::new(
                            "_unused".to_owned(),
                            summary_stats.iter(),
                        )
                        .compressed_with(utile::resource::Compression::Gzip)
                        .buffered()
                        .read()
                        .unwrap(),
                    )
                    .unwrap();
            })
            .join()
            .unwrap();
        }
    })
    .await
    .unwrap();

    log::info!("Cached summary stats for {}", phenotype.filename);

    drop(lock);
}

fn load_summary_stats(
    summary_stats: impl Iterator<Item = SummaryStats<hail::contig::GRCh37Contig>>,
    grch38: &mut FastaReader,
    grch37: &mut FastaReader,
    liftover: &liftover::LiftoverIndexed,
) -> Vec<SummaryStats<GRCh38Contig>> {
    log::info!("[load_summary_stats] Starting");

    let mut count = 0;

    let mut summary_stats: Vec<_> = summary_stats
        .inspect(|_| count += 1)
        .inspect_every(1_000_000, |i, _| {
            log::info!("[load_summary_stats] Processing summary stat {i}")
        })
        .filter(|s| {
            if s.beta_meta.is_none() {
                // log::warn!("Skipping summary stat: {:?}", s.at_range());
            }
            s.beta_meta.is_some()
        })
        .flat_map(move |s| {
            let at = s.at_range();

            {
                let found37 = grch37.query(&at).unwrap();

                if found37 != s.ref_allele {
                    log::warn!(
                        "[load_summary_stats] Reference mismatch: found37: {found37}, original: {}",
                        s.ref_allele
                    );
                }
            }

            liftover
                .map_range_raw(at.map_contig(|c| c.to_string()))
                .filter_map(|(mapped, was_flipped)| {
                    let mut new = s
                        .clone()
                        .map_contig(|_| GRCh38Contig::from_str(&mapped.name).unwrap());

                    // new.chr = mapped.name.clone();
                    new.pos = mapped.at.start + 1;
                    if was_flipped {
                        new.ref_allele = s.ref_allele.clone().reverse_complement(|b| b.complement())
                    }
                    if was_flipped {
                        new.alt = s.alt.clone().reverse_complement(|b| b.complement())
                    }

                    let found38 = match grch38.query(&mapped) {
                        Ok(found) => found,
                        Err(e) => {
                            log::warn!(
                                "[load_summary_stats] Error querying grch38: {e} ({mapped:?})"
                            );
                            return None;
                        }
                    };

                    match (found38 == new.ref_allele, found38 == new.alt) {
                        (true, true) => {
                            log::warn!("[load_summary_stats] Found both ref and alt in 38");
                            None
                        }
                        (true, false) => Some(new),
                        (false, true) => {
                            // log::warn!("Flipping ref/alt");
                            new.flip_ref_alt();
                            Some(new)
                        }
                        (false, false) => {
                            log::warn!("[load_summary_stats] Found neither ref nor alt in 38");
                            None
                        }
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect();

    log::info!("[load_summary_stats] {count} variants processed");

    summary_stats.sort_unstable_by(cmp_summary_stats);

    log::info!(
        "[load_summary_stats] {} variants sorted",
        summary_stats.len()
    );

    summary_stats
}

fn cmp_summary_stats(a: &SummaryStats<GRCh38Contig>, b: &SummaryStats<GRCh38Contig>) -> Ordering {
    Ord::cmp(
        &(a.at_range(), &a.ref_allele, &a.alt),
        &(b.at_range(), &b.ref_allele, &b.alt),
    )
}
