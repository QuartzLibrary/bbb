#![feature(map_try_insert)]
#![feature(let_chains)]
#![feature(iterator_try_collect)]
#![feature(mixed_integer_ops_unsigned_sub)]
#![feature(iter_map_windows)]

use std::{cmp::Ordering, collections::BTreeMap, sync::LazyLock};

use futures::{StreamExt, stream};
use genomes1000::simplified::SimplifiedRecord;
use itertools::Itertools;
use pan_ukbb::{PhenotypeManifestEntry, SummaryStats};
use tokio::task::spawn_blocking;
use utile::{
    cache::{FsCache, FsCacheEntry},
    plot::Histogram,
    resource::{RawResource, RawResourceExt},
};

use analysis::scores::Scores;

static CACHE: LazyLock<FsCache> = LazyLock::new(|| FsCache::new("/media/user/external2/data"));
// Used for intermediate files and avoid needing to keep all the variants in memory.
static TEMP_CACHE: LazyLock<FsCache> = LazyLock::new(|| CACHE.subfolder("temp"));

type FastaReader = biocore::fasta::IndexedFastaReader<std::io::BufReader<std::fs::File>>;

#[tokio::main]
async fn main() {
    const CHUNK_SIZE: usize = 10;
    const CONCURRENCY: usize = 1;

    env_logger::Builder::new()
        .filter_level(log::LevelFilter::Debug)
        .filter_module("reqwest", log::LevelFilter::Info)
        .filter_module("hyper_util", log::LevelFilter::Info)
        .filter_module("utile::resource", log::LevelFilter::Warn)
        .init();

    // let _clear_cache = ExecuteOnDrop::new(|| std::fs::remove_dir_all(&*TEMP_CACHE).unwrap());

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
            let exists = output_path(phenotype).as_ref().try_exists().unwrap();
            if exists {
                log::info!("[PanUKBB][{}] Found cached value", phenotype.filename);
            }
            !exists
        })
        .inspect(|phenotype| log::info!("[PanUKBB][{}] Processing", phenotype.filename))
        .chunks(CHUNK_SIZE);

    let phenotypes = phenotypes
        .into_iter()
        .map(|phenotypes| run_phenotypes(phenotypes.collect_vec(), lock, liftover))
        .map(|fut| async move { tokio::spawn(fut).await.unwrap() });

    stream::iter(phenotypes)
        .buffer_unordered(CONCURRENCY)
        .for_each(|()| async {})
        .await;

    log_memory(&mut sys);
}

#[tokio::test]
async fn compress() {
    fn output_path2(phenotype: &PhenotypeManifestEntry) -> FsCacheEntry {
        CACHE.entry(format!(
            "output/scores2/{}.json.bz",
            phenotype.filename.strip_suffix(".tsv.bgz").unwrap()
        ))
    }

    async fn run(phenotype: PhenotypeManifestEntry) {
        use utile::resource::RawResource;

        let entry = output_path(&phenotype);
        let entry2 = output_path2(&phenotype);
        if !entry2.try_exists().unwrap() {
            log::info!("[PanUKBB][{}] Compressing", phenotype.filename);
            entry2
                .write_file_async(
                    entry
                        .compressed_with(utile::resource::Compression::Gzip)
                        .buffered()
                        .read_async()
                        .await
                        .unwrap(),
                )
                .await
                .unwrap();
        }
    }
    env_logger::Builder::new()
        .filter_level(log::LevelFilter::Debug)
        .filter_module("reqwest", log::LevelFilter::Info)
        .filter_module("hyper_util", log::LevelFilter::Info)
        .filter_module("utile::resource", log::LevelFilter::Warn)
        .init();

    let mut sys = sysinfo::System::new_all();

    log_memory(&mut sys);

    let phenotypes = PhenotypeManifestEntry::load_default()
        .await
        .unwrap()
        .into_iter()
        .filter(analysis::util::passes_qc)
        .filter(|phenotype| {
            let exists = output_path(phenotype).try_exists().unwrap();
            if !exists {
                log::info!("[PanUKBB][{}] Missing value", phenotype.filename);
            }
            exists
        })
        .inspect(|phenotype| log::info!("[PanUKBB][{}] Processing", phenotype.filename))
        .map(|phenotype| run(phenotype))
        .map(|fut| async move { tokio::spawn(fut).await.unwrap() });

    stream::iter(phenotypes).for_each(|a| a).await;

    log_memory(&mut sys);
}

async fn run_phenotypes(
    phenotypes: Vec<PhenotypeManifestEntry>,
    lock: &'static tokio::sync::Mutex<()>,
    liftover: &'static liftover::LiftoverIndexed,
) {
    for p in &phenotypes {
        let path = output_path(p);
        std::fs::create_dir_all(path.as_ref().parent().unwrap()).unwrap();
    }

    log::info!("[PanUKBB][{}] Loading summary stats", phenotypes.len());

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

    let all_summary_stats = itertools::kmerge_by(
        all_summary_stats,
        |a: &(String, SummaryStats), b: &(String, SummaryStats)| {
            let (_, a) = a;
            let (_, b) = b;
            cmp_summary_stats(a, b) == Ordering::Less
        },
    );

    log::info!("[PanUKBB][{}] Scoring", phenotypes.len());

    let mut multi_scores = scores(&phenotypes, all_summary_stats).await;

    log::info!("[PanUKBB][{}] Done scoring", phenotypes.len());

    // for phenotype in &phenotypes {
    //     summary_stats_temp_path(phenotype)
    //         .invalidate_async()
    //         .await
    //         .unwrap();
    // }

    for scores in multi_scores.raw_scores() {
        Histogram {
            data: scores.scores().map(|s| s.score).collect(),
            bins: Some(100),
        }
        .show();
    }

    multi_scores.finalize();

    for score in multi_scores.all() {
        log::info!("[PanUKBB][{}] Writing", score.phenotype.filename);

        let path = output_path(&score.phenotype);
        path.write_json(score).unwrap();

        // output_path(&score.phenotype).write_file_with(|file| {
        //     let file = brotli::CompressorWriter::new(file, 4096, 9, 20);
        //     serde_json::to_writer(file, score).unwrap();
        //     Ok(())
        // })
        // .unwrap();

        log::info!("[PanUKBB][{}] Finished scoring", score.phenotype.filename);
    }
}

async fn scores(
    phenotypes: &[PhenotypeManifestEntry],
    summary_stats: impl Iterator<Item = (String, SummaryStats)>,
) -> MultiScores {
    let mut summary_stats = summary_stats.peekable();
    let mut scores = MultiScores::new(phenotypes).await;

    let (sample_names, variants) = genomes1000::load_all_simplified().await;
    let mut variants = variants.peekable();

    while let Some((_filename, summary_stat)) = summary_stats.peek() {
        while let Some(record) = variants.peek() {
            let cmp = cmp_variant(record, summary_stat);
            match cmp {
                Ordering::Equal => {
                    let (filename, summary_stat) = summary_stats.next().unwrap();

                    scores.push_variant(&filename, summary_stat, &sample_names, record);

                    // variants.next(); // We don't consume the variant until all summary stats have been processed!
                    break;
                }
                Ordering::Less => {
                    variants.next();
                }
                Ordering::Greater => {
                    let (filename, summary_stat) = summary_stats.next().unwrap();

                    scores.push_missing(&filename, summary_stat);

                    break;
                }
            }
        }
    }

    scores
}

fn cmp_variant(a: &SimplifiedRecord, b: &SummaryStats) -> Ordering {
    Ord::cmp(&a.at(), &b.at())
        .then_with(|| Ord::cmp(&a.reference_allele, &b.ref_allele))
        .then_with(|| Ord::cmp(&a.alternate_allele, &b.alt))
}

async fn summary_stats(
    phenotype: PhenotypeManifestEntry,
    liftover: &'static liftover::LiftoverIndexed,
    lock: &'static tokio::sync::Mutex<()>,
) -> impl Iterator<Item = SummaryStats> {
    let summary_stats_entry = summary_stats_temp_path(&phenotype);

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

    let summary_stats = SummaryStats::load(
        phenotype
            .summary_stats_resource()
            .log_progress()
            .with_fs_cache(&CACHE)
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
    summary_stats: impl Iterator<Item = SummaryStats>,
    grch38: &mut FastaReader,
    grch37: &mut FastaReader,
    liftover: &liftover::LiftoverIndexed,
) -> Vec<SummaryStats> {
    log::info!("[load_summary_stats] Starting");

    let mut count = 0;

    let mut summary_stats: Vec<_> = summary_stats
        .enumerate()
        .inspect(|(i, _)| {
            count += 1;
            if i % 1_000_000 == 0 {
                log::info!("Processing variant {i}");
            }
        })
        .map(|(_, v)| v)
        .filter(|s| {
            if s.beta_meta.is_none() {
                // log::warn!("Skipping variant: {:?}", s.at_range());
            }
            s.beta_meta.is_some()
        })
        .flat_map(move |s| {
            let at = s.at_range();

            {
                let found37 = grch37.query(&at).unwrap();

                if found37 != s.ref_allele {
                    log::warn!(
                        "Reference mismatch: found37: {found37}, original: {}",
                        s.ref_allele
                    );
                }
            }

            liftover
                .map_range_raw(at)
                .filter_map(|(mapped, was_flipped)| {
                    let mut new = s.clone();

                    new.chr = mapped.name.clone();
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
                            log::warn!("Error querying grch38: {e} ({mapped:?})");
                            return None;
                        }
                    };

                    match (found38 == new.ref_allele, found38 == new.alt) {
                        (true, true) => {
                            log::warn!("Found both ref and alt in 38");
                            None
                        }
                        (true, false) => Some(new),
                        (false, true) => {
                            // log::warn!("Flipping ref/alt");
                            new.flip_ref_alt();
                            Some(new)
                        }
                        (false, false) => {
                            log::warn!("Found neither ref nor alt in 38");
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

fn cmp_summary_stats(a: &SummaryStats, b: &SummaryStats) -> Ordering {
    Ord::cmp(
        &(a.at_range(), &a.ref_allele, &a.alt),
        &(b.at_range(), &b.ref_allele, &b.alt),
    )
}

async fn load_liftover() -> liftover::LiftoverIndexed {
    spawn_blocking(|| {
        log::info!("[Liftover] Loading");
        let liftover = liftover::Liftover::load(
            hail::resource::HailCommonResource::grch37_to_grch38_liftover_chain()
                .log_progress()
                .with_global_fs_cache(),
        )
        .unwrap()
        .indexed();
        log::info!("[Liftover] Loaded");
        liftover
    })
    .await
    .unwrap()
}

fn summary_stats_temp_path(phenotype: &PhenotypeManifestEntry) -> FsCacheEntry {
    TEMP_CACHE.entry(phenotype.filename.clone())
}

fn output_path(phenotype: &PhenotypeManifestEntry) -> FsCacheEntry {
    CACHE.entry(format!("output/scores/{}.json", phenotype.filename))
}

fn log_memory(sys: &mut sysinfo::System) {
    fn gb(bytes: u64) -> f64 {
        bytes as f64 / 1024.0 / 1024.0 / 1024.0
    }
    log::info!("=> system:");
    // RAM and swap information:
    log::info!("total memory: {}GB", gb(sys.total_memory()));
    log::info!("used memory : {}GB", gb(sys.used_memory()));
    log::info!("total swap  : {}GB", gb(sys.total_swap()));
    log::info!("used swap   : {}GB", gb(sys.used_swap()));
    log::info!(
        "process memory: {}GB",
        gb(sys
            .process(sysinfo::Pid::from_u32(std::process::id()))
            .unwrap()
            .memory())
    );
}

pub struct MultiScores {
    pub all: BTreeMap<String, Scores>,
}

impl MultiScores {
    pub async fn new(phenotypes: &[PhenotypeManifestEntry]) -> Self {
        let mut all = BTreeMap::new();
        for p in phenotypes {
            all.insert(p.filename.clone(), Scores::new(p.clone()).await);
        }
        Self { all }
    }
    pub fn push_missing(&mut self, filename: &str, summary_stat: SummaryStats) {
        self.all
            .get_mut(filename)
            .unwrap()
            .push_missing(summary_stat);
    }
    pub fn push_variant(
        &mut self,
        filename: &str,
        summary_stat: SummaryStats,
        sample_names: &[String],
        record: &SimplifiedRecord,
    ) {
        self.all
            .get_mut(filename)
            .unwrap()
            .push_variant(summary_stat, sample_names, record);
    }

    pub fn finalize(&mut self) {
        for scores in self.all.values_mut() {
            scores.finalize();
        }
    }

    pub fn all(&self) -> impl Iterator<Item = &Scores> {
        self.all.values()
    }

    pub fn raw_scores(&self) -> impl Iterator<Item = &Scores> {
        self.all.values()
    }
}
