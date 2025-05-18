#![feature(let_chains)]

mod shared;

use std::{cmp::Ordering, collections::BTreeSet};

use futures::{StreamExt, stream};
use genomes1000::{Genomes1000Fs, simplified::SimplifiedRecord};
use hail::contig::GRCh38Contig;
use ordered_float::NotNan;
use pan_ukbb::{PhenotypeManifestEntry, SummaryStats};
use shared::{output_path_pvalues, summary_stats_pvalues_temp_path};
use utile::resource::{RawResource, RawResourceExt};

use analysis::{
    pvalues::PhenotypeTopPValues,
    util::{FirstN, InspectEvery},
};

const TOP_PVALUES_N: usize = 3_000;

use self::shared::{cmp_variant, load_liftover, log_memory, summary_stats::summary_stats};

#[tokio::main]
async fn main() {
    const CONCURRENCY: usize = 50;

    env_logger::Builder::new()
        .filter_level(log::LevelFilter::Debug)
        .filter_module("reqwest", log::LevelFilter::Info)
        .filter_module("hyper_util", log::LevelFilter::Info)
        .filter_module("utile::resource", log::LevelFilter::Warn)
        .init();

    // let _clear_cache =
    //     ExecuteOnDrop::new(|| std::fs::remove_dir_all(&*shared::TEMP_CACHE).unwrap());

    let mut sys = sysinfo::System::new_all();

    log_memory(&mut sys);

    let lock = &*Box::leak(Box::new(tokio::sync::Mutex::new(())));
    let liftover = &*Box::leak(Box::new(load_liftover().await));

    let phenotypes = PhenotypeManifestEntry::load_default()
        .await
        .unwrap()
        .into_iter()
        .filter(analysis::util::passes_qc)
        .filter(|phenotype| {
            let exists = output_path_pvalues(phenotype).try_exists().unwrap();
            if exists {
                log::info!(
                    "[PanUKBB][pvalues][{}] Found cached value",
                    phenotype.filename
                );
            }
            !exists
        })
        .inspect(|phenotype| log::info!("[PanUKBB][pvalues][{}] Processing", phenotype.filename))
        .map(|p| run_pvalues(p, lock, liftover))
        .map(|fut| async move { tokio::spawn(fut).await.unwrap() });

    stream::iter(phenotypes)
        .buffer_unordered(CONCURRENCY)
        .for_each(|()| async {})
        .await;

    log_memory(&mut sys);
}

async fn run_pvalues(
    phenotype: PhenotypeManifestEntry,
    lock: &'static tokio::sync::Mutex<()>,
    liftover: &'static liftover::LiftoverIndexed,
) {
    let path = output_path_pvalues(&phenotype);
    std::fs::create_dir_all(path.as_ref().parent().unwrap()).unwrap();

    log::info!(
        "[PanUKBB][pvalues][{}] Preparing top pvalues",
        phenotype.filename
    );

    let top_pvalues = load_and_cache_top_pvalues(phenotype.clone(), lock, liftover).await;

    let matched_stats = match_summary_stats(top_pvalues)
        .await
        .inspect_every(1_000, |i, _| {
            log::info!("[pvalues] Processing matched summary stat {i}")
        });

    let (sample_names, _) = genomes1000::load_contig(genomes1000::GRCh38Contig::CHR22)
        .await
        .unwrap();

    let mut result = PhenotypeTopPValues::new(phenotype.clone()).await;
    for (stat, record) in matched_stats {
        match record {
            Some(record) => {
                result.push_variant(stat, &sample_names, &record);
            }
            None => {
                result.push_missing(stat);
            }
        }
    }

    log::info!("[PanUKBB][pvalues][{}] Writing", phenotype.filename);

    output_path_pvalues(&phenotype)
        .write_file_with(|file| {
            let file = brotli::CompressorWriter::new(file, 4096, 9, 20);
            serde_json::to_writer(file, &result).unwrap();
            Ok(())
        })
        .unwrap();

    log::info!(
        "[PanUKBB][pvalues][{}] Finished scoring",
        phenotype.filename
    );
}

async fn match_summary_stats(
    summary_stats: impl Iterator<Item = SummaryStats<GRCh38Contig>>,
) -> impl Iterator<Item = (SummaryStats<GRCh38Contig>, Option<SimplifiedRecord>)> {
    let mut genomes1000 = Genomes1000Fs::new().await.unwrap();

    summary_stats.map(move |summary_stat| {
        for record in genomes1000
            .query_simplified(
                &summary_stat
                    .at_range()
                    .map_contig(|c| genomes1000::GRCh38Contig::new(c.as_ref()).unwrap()),
            )
            .unwrap()
        {
            match cmp_variant(&record, &summary_stat) {
                Ordering::Equal => {
                    return (summary_stat, Some(record));
                }
                Ordering::Less => {}
                Ordering::Greater => {
                    return (summary_stat, None);
                }
            }
        }

        (summary_stat, None)
    })
}

async fn load_and_cache_top_pvalues(
    phenotype: PhenotypeManifestEntry,
    lock: &'static tokio::sync::Mutex<()>,
    liftover: &'static liftover::LiftoverIndexed,
) -> impl Iterator<Item = SummaryStats<GRCh38Contig>> {
    let path = summary_stats_pvalues_temp_path(&phenotype);

    let exists = path.try_exists().unwrap();

    if !exists {
        let mut max_pvalue = None;
        let mut max_pvalue_hq = None;

        let mut min_pvalue = None;
        let mut min_pvalue_hq = None;

        let mut all_pvalues = vec![];
        let mut all_pvalues_hq = vec![];

        let mut pvalues = FirstN::new(TOP_PVALUES_N, |a: &SummaryStats<GRCh38Contig>, b| {
            Ord::cmp(
                &a.neglog10_pval_meta.unwrap(),
                &b.neglog10_pval_meta.unwrap(),
            )
            .reverse()
        });
        let mut pvalues_hq = FirstN::new(TOP_PVALUES_N, |a: &SummaryStats<GRCh38Contig>, b| {
            Ord::cmp(
                &a.neglog10_pval_meta_hq.unwrap(),
                &b.neglog10_pval_meta_hq.unwrap(),
            )
            .reverse()
        });

        for stat in summary_stats(phenotype.clone(), liftover, lock)
            .await
            .inspect_every(1_000_000, |i, _| {
                log::info!("[pvalues] Processing summary stat {i}")
            })
        {
            if let Some(neglog10_pval_meta) = stat.neglog10_pval_meta
                && stat.beta_meta.is_some()
            {
                all_pvalues.push(neglog10_pval_meta);
                max_pvalue = Some(
                    max_pvalue
                        .unwrap_or(neglog10_pval_meta)
                        .max(neglog10_pval_meta),
                );
                min_pvalue = Some(
                    min_pvalue
                        .unwrap_or(neglog10_pval_meta)
                        .min(neglog10_pval_meta),
                );
                pvalues.push_ref(&stat);
            }
            if let Some(neglog10_pval_meta_hq) = stat.neglog10_pval_meta_hq
                && stat.beta_meta_hq.is_some()
            {
                all_pvalues_hq.push(neglog10_pval_meta_hq);
                max_pvalue_hq = Some(
                    max_pvalue_hq
                        .unwrap_or(neglog10_pval_meta_hq)
                        .max(neglog10_pval_meta_hq),
                );
                min_pvalue_hq = Some(
                    min_pvalue_hq
                        .unwrap_or(neglog10_pval_meta_hq)
                        .min(neglog10_pval_meta_hq),
                );
                pvalues_hq.push_ref(&stat);
            }
        }

        // Histogram {
        //     data: all_pvalues,
        //     bins: Some(100),
        // }
        // .show();
        // Histogram {
        //     data: all_pvalues_hq,
        //     bins: Some(100),
        // }
        // .show();

        log::info!(
            "[PanUKBB][pvalues][{}] Max pvalue: {max_pvalue:?}, Min pvalue: {min_pvalue:?}, count: {}",
            phenotype.filename,
            pvalues.len(),
        );
        log::info!(
            "[PanUKBB][pvalues][{}] Max pvalue HQ: {max_pvalue_hq:?}, Min pvalue HQ: {min_pvalue_hq:?}, count: {}",
            phenotype.filename,
            pvalues_hq.len(),
        );

        log::info!(
            "[PanUKBB][pvalues][{}] Matching summary stats",
            phenotype.filename
        );

        assert_eq!(
            pvalues.first().unwrap().neglog10_pval_meta.unwrap(),
            max_pvalue.unwrap()
        );
        assert_eq!(
            pvalues_hq.first().unwrap().neglog10_pval_meta_hq.unwrap(),
            max_pvalue_hq.unwrap()
        );

        let all_stats: BTreeSet<_> = pvalues.into_iter().chain(pvalues_hq).collect();

        log::info!(
            "[PanUKBB][pvalues][{}] Found {} summary stats",
            phenotype.filename,
            all_stats.len()
        );

        path.write_file(
            utile::resource::iter::IterToJsonLinesResource::new(
                "_unused".to_owned(),
                all_stats.iter(),
            )
            .compressed_with(utile::resource::Compression::Gzip)
            .buffered()
            .read()
            .unwrap(),
        )
        .unwrap();
    }

    let mut all: Vec<SummaryStats<GRCh38Contig>> = path
        .decompressed_with(utile::resource::Compression::Gzip)
        .read_json_lines()
        .unwrap()
        .map(|s| s.unwrap())
        .collect();

    // for s in &all {
    //     assert!(s.beta_meta.is_some());
    //     assert!(s.beta_meta_hq.is_some());
    // }

    let neg_infinity = NotNan::new(f64::NEG_INFINITY).unwrap();

    let pvalues = {
        all.sort_unstable_by_key(|s| -(s.neglog10_pval_meta.unwrap_or(neg_infinity)));
        all[..TOP_PVALUES_N].to_vec()
    };
    log::info!("Pvalues: {:?}", pvalues.len());
    // log::warn!("First pvalue: {:?}", &pvalues[0]);
    // log::warn!("Last pvalue: {:?}", &pvalues[TOP_PVALUES_N - 1]);

    let pvalues_hq = {
        all.sort_unstable_by_key(|s| -(s.neglog10_pval_meta_hq.unwrap_or(neg_infinity)));
        all[..TOP_PVALUES_N].to_vec()
    };
    log::info!("Pvalues HQ: {:?}", pvalues_hq.len());
    // log::warn!("First pvalue HQ: {:?}", &pvalues_hq[0]);
    // log::warn!("Last pvalue HQ: {:?}", &pvalues_hq[pvalues_hq.len() - 1]);

    let done: BTreeSet<_> = pvalues.into_iter().chain(pvalues_hq).collect();

    if done.len() != all.len() {
        log::info!("Done: {:?}", done.len());
        log::info!("All: {:?}", all.len());
        // for s in &all {
        //     if !done.contains(&s) {
        //         // println!("{:?}", );
        //         log::warn!("Missing: {:?}", &s);
        //         assert_eq!(done.len(), all.len());
        //     }
        // }
        // assert_eq!(done.len(), all.len());
    }

    done.into_iter()
}
