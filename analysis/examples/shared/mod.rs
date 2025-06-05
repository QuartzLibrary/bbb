#![allow(dead_code)] // There is code used in only one of the two binaries.

pub mod summary_stats;

use std::{cmp::Ordering, str::FromStr, sync::LazyLock};

use genomes1000::simplified::SimplifiedRecord;
use hail::contig::GRCh38Contig;
use pan_ukbb::{PhenotypeManifestEntry, SummaryStats};
use tokio::task::spawn_blocking;
use utile::{
    cache::{FsCache, FsCacheEntry},
    resource::RawResourceExt,
};

/// Where we put the putput
pub static DATA_FOLDER: LazyLock<FsCache> =
    LazyLock::new(|| FsCache::new("/media/user/external2/data"));
pub static OUTPUT_FOLDER: LazyLock<FsCache> = LazyLock::new(|| DATA_FOLDER.subfolder("output"));

// Used for intermediate files and avoid needing to keep all the variants in memory.
pub static TEMP_CACHE: LazyLock<FsCache> = LazyLock::new(|| DATA_FOLDER.subfolder("temp"));

pub fn output_path(phenotype: &PhenotypeManifestEntry) -> FsCacheEntry {
    let name = phenotype.filename.strip_suffix(".tsv.bgz").unwrap();
    OUTPUT_FOLDER.entry(format!("scores/{name}.json.br"))
}
pub fn output_path_pvalues(phenotype: &PhenotypeManifestEntry) -> FsCacheEntry {
    let name = phenotype.filename.strip_suffix(".tsv.bgz").unwrap();
    OUTPUT_FOLDER.entry(format!("pvalues/{name}.json.br"))
}

/// Intermediate file for the lifted-over and re-sorted summary stats.
pub fn summary_stats_temp_path(phenotype: &PhenotypeManifestEntry) -> FsCacheEntry {
    let filename = phenotype.filename.clone();
    TEMP_CACHE.entry(format!("lifted/{filename}"))
}
/// Intermediate file for the top summary stats by pvalue.
pub fn summary_stats_pvalues_temp_path(phenotype: &PhenotypeManifestEntry) -> FsCacheEntry {
    let filename = phenotype.filename.clone();
    TEMP_CACHE.entry(format!("pvalues/{filename}"))
}

pub fn cmp_variant(a: &SimplifiedRecord, b: &SummaryStats<GRCh38Contig>) -> Ordering {
    Ord::cmp(
        &a.at(),
        &b.at()
            .map_contig(|c| genomes1000::GRCh38Contig::from_str(c.as_ref()).unwrap()),
    )
    .then_with(|| Ord::cmp(&a.reference_allele, &b.ref_allele))
    .then_with(|| Ord::cmp(&a.alternate_allele, &b.alt))
}

pub async fn load_liftover() -> liftover::LiftoverIndexed {
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

pub fn log_memory(sys: &mut sysinfo::System) {
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
