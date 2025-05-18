mod shared;

use std::{cmp::Ordering, collections::BTreeMap};

use futures::{StreamExt, stream};
use genomes1000::simplified::SimplifiedRecord;
use hail::contig::GRCh38Contig;
use itertools::Itertools;
use pan_ukbb::{PhenotypeManifestEntry, SummaryStats};
use utile::plot::Histogram;

use analysis::{scores::Scores, util::InspectEvery};

use self::shared::{
    cmp_variant, load_liftover, log_memory, output_path, summary_stats::merged_summary_stats,
};

#[tokio::main]
async fn main() {
    const CHUNK_SIZE: usize = 1;
    const CONCURRENCY: usize = 15;

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
            let exists = output_path(phenotype).try_exists().unwrap();
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

    let all_summary_stats = merged_summary_stats(phenotypes.clone(), lock, liftover).await;

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

        output_path(&score.phenotype)
            .write_file_with(|file| {
                let file = brotli::CompressorWriter::new(file, 4096, 9, 20);
                serde_json::to_writer(file, score).unwrap();
                Ok(())
            })
            .unwrap();

        log::info!("[PanUKBB][{}] Finished scoring", score.phenotype.filename);
    }
}

async fn scores(
    phenotypes: &[PhenotypeManifestEntry],
    summary_stats: impl Iterator<Item = (String, SummaryStats<GRCh38Contig>)>,
) -> MultiScores {
    let mut summary_stats = summary_stats.peekable();
    let mut scores = MultiScores::new(phenotypes).await;

    let (sample_names, variants) = genomes1000::load_all_simplified().await;
    let mut variants = variants
        .inspect_every(1_000_000, |i, _| {
            log::info!("[scores] Processing variant {i}")
        })
        .peekable();

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
    pub fn push_missing(&mut self, filename: &str, summary_stat: SummaryStats<GRCh38Contig>) {
        self.all
            .get_mut(filename)
            .unwrap()
            .push_missing(summary_stat);
    }
    pub fn push_variant(
        &mut self,
        filename: &str,
        summary_stat: SummaryStats<GRCh38Contig>,
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
