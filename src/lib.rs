#![feature(btree_set_entry)]
#![feature(let_chains)]

pub mod send_option;
pub mod util;

use futures::{
    stream::{AbortHandle, Abortable},
    StreamExt,
};
use gloo_worker::{HandlerId, Worker, WorkerScope};
use ordered_float::NotNan;
use serde::{de::DeserializeOwned, Deserialize, Serialize};
use std::{
    cell::RefCell,
    cmp::Ordering,
    collections::BTreeMap,
    fmt::Debug,
    future::Future,
    io,
    ops::{Deref, DerefMut, Range},
    sync::OnceLock,
    thread::LocalKey,
    time::Duration,
};

use hail::contig::GRCh38Contig;
use pan_ukbb::{PhenotypeManifestEntry, SummaryStats};
use utile::{
    collections::counting_set::CountingBTreeSet,
    resource::{Compression, RawResource, RawResourceExt, UrlResource},
};

use analysis::{
    pvalues::{PhenotypeTopPValues, SampleGenotype, TopPValueVariant},
    scores::{Scores, Stats},
    util::SummaryStatKey,
};

thread_local! {
    static SCORES_TASK: RefCell<Option<(HandlerIdOrd, AbortHandle)>> = const { RefCell::new(None) };
    static EDIT_ANALYSIS_TASK: RefCell<Option<(HandlerIdOrd, AbortHandle)>> = const { RefCell::new(None) };
}

static SCORES: OnceLock<Scores> = OnceLock::new();
static PVALUES: OnceLock<PhenotypeTopPValues> = OnceLock::new();

async fn scores(id: HandlerId) -> &'static Scores {
    wait_for(id, || SCORES.get()).await
}
async fn pvalues(id: HandlerId) -> &'static PhenotypeTopPValues {
    wait_for(id, || PVALUES.get()).await
}

async fn stats(id: HandlerId, use_hq: bool) -> Option<Stats> {
    static STATS: OnceLock<Option<Stats>> = OnceLock::new();
    static STATS_HQ: OnceLock<Option<Stats>> = OnceLock::new();

    let store = if use_hq { &STATS_HQ } else { &STATS };

    if let Some(s) = store.get() {
        return *s;
    }

    let scores = scores(id).await;
    let s = if use_hq {
        scores.stats_hq()
    } else {
        scores.stats()
    };
    store.set(s).unwrap();
    s
}
/// The variants with the highest p-values.
async fn top_pvalue_variants(id: HandlerId, use_hq: bool) -> &'static [&'static TopPValueVariant] {
    fn get_top_variants_hq(pvalues: &PhenotypeTopPValues) -> Vec<&TopPValueVariant> {
        let mut top_variants: Vec<_> = pvalues
            .top_variants
            .values()
            .filter(|v| v.stat.neglog10_pval_meta_hq.is_some() && v.stat.beta_meta_hq.is_some())
            .collect();
        top_variants.sort_unstable_by_key(|v| -v.stat.neglog10_pval_meta_hq.unwrap());
        top_variants
    }
    fn get_top_variants(pvalues: &PhenotypeTopPValues) -> Vec<&TopPValueVariant> {
        let mut top_variants: Vec<_> = pvalues
            .top_variants
            .values()
            .filter(|v| v.stat.neglog10_pval_meta.is_some() && v.stat.beta_meta.is_some())
            .collect();
        top_variants.sort_unstable_by_key(|v| -v.stat.neglog10_pval_meta.unwrap());
        top_variants
    }

    static TOP: OnceLock<Vec<&'static TopPValueVariant>> = OnceLock::new();
    static TOP_HQ: OnceLock<Vec<&'static TopPValueVariant>> = OnceLock::new();

    let store = if use_hq { &TOP_HQ } else { &TOP };

    if let Some(top) = store.get() {
        return top.as_slice();
    }

    log::info!("[Worker][{id:?}] Getting top p-values (use_hq: {use_hq}).");
    let start = web_time::Instant::now();

    let pvalues = pvalues(id).await;
    let top = if use_hq {
        get_top_variants_hq(pvalues)
    } else {
        get_top_variants(pvalues)
    };
    store.set(top).unwrap();

    log::info!(
        "[Worker][{id:?}] Got top p-values in {:?}ms",
        start.elapsed().as_millis()
    );

    store.get().unwrap().as_slice()
}

pub struct WorkerStruct {}

#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
pub enum Input {
    Set {
        origin: String,
        file: String,
    },
    GetScores {
        edit_count: isize,
        normalise: bool,
        use_hq: bool,
        top_pvalues: usize,
    },
    GetEditAnalysis {
        edit_count: isize,
        normalise: bool,
        use_hq: bool,
        top_pvalues: usize,
    },
}

#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
#[expect(clippy::large_enum_variant)]
pub enum Output {
    Init(OutputInit),
    GetScores(Option<OutputGetScores>),
    GetEditAnalysis(Option<OutputGetEditAnalysis>),
}
#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
pub struct OutputInit {
    // Horrible hack (AsJson)
    pub phenotype: AsJson<PhenotypeManifestEntry>,

    pub sample_count: usize,

    pub stats: Option<Stats>,
    pub stats_hq: Option<Stats>,

    pub full_score_range: Option<Range<f64>>,
    pub full_score_range_hq: Option<Range<f64>>,

    pub norm_stats: Option<Stats>,
    pub norm_stats_hq: Option<Stats>,

    pub norm_full_score_range: Option<Range<f64>>,
    pub norm_full_score_range_hq: Option<Range<f64>>,
}
#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
pub struct OutputGetScores {
    pub scores: Vec<f64>,

    pub stats: Stats,
}
impl OutputGetScores {
    pub fn normalised(self, mean: f64, std_dev: f64) -> Self {
        Self {
            scores: self
                .scores
                .into_iter()
                .map(|v| (v - mean) / std_dev)
                .collect(),
            stats: self.stats.normalised(mean, std_dev),
        }
    }
}

#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
pub struct OutputGetEditAnalysis {
    pub use_hq: bool,
    /// Pre-sorted by p-value.
    pub edits: Vec<VariantInfo>,
}
impl OutputGetEditAnalysis {
    pub fn normalised(self, std_dev: NotNan<f64>, std_dev_hq: NotNan<f64>) -> Self {
        Self {
            use_hq: self.use_hq,
            edits: self
                .edits
                .into_iter()
                .map(|v| v.normalised(std_dev, std_dev_hq))
                .collect(),
        }
    }
}
#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
pub struct VariantInfo {
    pub count: usize,
    pub stat: AsJson<SummaryStats<GRCh38Contig>>,
    pub ploidy_dosage: CountingBTreeSet<VariantSampleInfo>,
}
impl VariantInfo {
    pub fn new(key: SummaryStatKey, scores: &Scores) -> Self {
        let SummaryStatKey {
            chr,
            pos,
            ref_allele,
            alt,
            ..
        } = key;

        let stat = scores
            .summary_stats
            .iter()
            .find(|s| s.chr == chr && s.pos == pos && s.ref_allele == ref_allele && s.alt == alt)
            .unwrap()
            .clone();

        VariantInfo {
            count: 0,
            stat: AsJson(stat),
            ploidy_dosage: CountingBTreeSet::new(),
        }
    }
    pub fn normalised(self, std_dev: NotNan<f64>, std_dev_hq: NotNan<f64>) -> Self {
        Self {
            count: self.count,
            stat: AsJson(self.stat.0.normalised(std_dev, std_dev_hq)),
            ploidy_dosage: self.ploidy_dosage,
        }
    }
}
#[derive(Debug, Clone, Eq, PartialEq, PartialOrd, Ord, Hash)]
#[derive(Serialize, Deserialize)]
pub struct VariantSampleInfo {
    pub dosage: u8,
    pub ploidy: u8,
}

impl Worker for WorkerStruct {
    type Input = Input;
    type Message = ();
    type Output = Output;

    fn create(_scope: &WorkerScope<Self>) -> Self {
        Self {}
    }

    fn update(&mut self, _scope: &WorkerScope<Self>, _msg: Self::Message) {
        unreachable!()
    }

    fn received(&mut self, scope: &WorkerScope<Self>, msg: Self::Input, id: HandlerId) {
        log::info!("[Worker][{id:?}] Received: {msg:?}");
        match msg {
            Input::Set { file, origin } => {
                let scope = scope.clone();
                let respond = move |scores: &Scores| {
                    let stats = scores.stats();
                    let stats_hq = scores.stats_hq();

                    let full_score_range = scores.full_score_range();
                    let full_score_range_hq = scores.full_score_range_hq();

                    let norm_stats = stats.map(|s| s.normalised_self());
                    let norm_stats_hq = stats_hq.map(|s| s.normalised_self());

                    let norm_full_score_range = stats.map(|s| {
                        let full_score_range = full_score_range.clone().unwrap();
                        s.normalise_value(full_score_range.start)
                            ..s.normalise_value(full_score_range.end)
                    });
                    let norm_full_score_range_hq = stats_hq.map(|s| {
                        let full_score_range_hq = full_score_range_hq.clone().unwrap();
                        s.normalise_value(full_score_range_hq.start)
                            ..s.normalise_value(full_score_range_hq.end)
                    });

                    let response = Output::Init(OutputInit {
                        phenotype: AsJson(scores.phenotype.clone()),

                        sample_count: scores.scores.len(),

                        stats,
                        stats_hq,

                        full_score_range,
                        full_score_range_hq,

                        norm_stats,
                        norm_stats_hq,

                        norm_full_score_range,
                        norm_full_score_range_hq,
                    });

                    scope.respond(id, response);
                };
                wasm_bindgen_futures::spawn_local(async move {
                    let scores = fetch_scores(id, &origin, file.clone()).await.unwrap();
                    SCORES.set(scores.clone()).unwrap();

                    respond(SCORES.get().unwrap());

                    // Give it a second for the initial scores to be sent.
                    utile::time::sleep(Duration::from_millis(1000)).await;

                    let pvalues = fetch_pvalues(id, &origin, file.clone()).await.unwrap();
                    PVALUES.set(pvalues).unwrap();
                });
            }
            Input::GetScores {
                edit_count,
                normalise,
                use_hq,
                top_pvalues,
            } => {
                let scope = scope.clone();
                let f = async move {
                    let mut data = compute_scores(id, edit_count, use_hq, top_pvalues).await;
                    if normalise && let Some(s) = stats(id, use_hq).await {
                        data = data.map(|d| d.normalised(s.mean, s.std_dev));
                    }
                    log::info!("[Worker][{id:?}] Sending scores.");
                    scope.respond(id, Output::GetScores(data))
                };
                spawn_task(&SCORES_TASK, id, f);
            }
            Input::GetEditAnalysis {
                edit_count,
                normalise,
                use_hq,
                top_pvalues,
            } => {
                let scope = scope.clone();
                let f = async move {
                    let mut data = compute_edit_analysis(id, edit_count, use_hq, top_pvalues).await;
                    if let Some(data) = &mut data {
                        data.edits.sort_unstable_by_key(|v| {
                            if use_hq {
                                -v.stat.neglog10_pval_meta_hq.unwrap()
                            } else {
                                -v.stat.neglog10_pval_meta.unwrap()
                            }
                        });
                    }
                    if normalise
                        && let Some(s) = stats(id, false).await
                        && let Some(s_hq) = stats(id, true).await
                    {
                        let std_dev = NotNan::new(s.std_dev).unwrap();
                        let std_dev_hq = NotNan::new(s_hq.std_dev).unwrap();
                        data = data.map(|d| d.normalised(std_dev, std_dev_hq));
                    }
                    log::info!("[Worker][{id:?}] Sending edit analysis.");
                    scope.respond(id, Output::GetEditAnalysis(data))
                };
                spawn_task(&EDIT_ANALYSIS_TASK, id, f);
            }
        }
    }
}

/// `async` to give it a chance to be interrupted.
async fn compute_scores(
    id: HandlerId,
    edit_count: isize,
    use_hq: bool,
    top_pvalues: usize,
) -> Option<OutputGetScores> {
    log::info!("[Worker][{id:?}] Computing scores.");

    let values = compute_scores_(id, edit_count, use_hq, top_pvalues).await?;

    let iter = || values.iter().copied();
    let mean = analysis::util::mean(iter()).unwrap();
    let std_dev = analysis::util::std_dev(mean, iter()).unwrap();
    let iter = || values.iter().copied().map(|v| NotNan::new(v).unwrap());
    let min = *iter().min().unwrap();
    let max = *iter().max().unwrap();

    Some(OutputGetScores {
        scores: values,
        stats: Stats {
            mean,
            std_dev,
            min,
            max,
        },
    })
}
/// `async` to give it a chance to be interrupted.
async fn compute_scores_(
    id: HandlerId,
    edit_count: isize,
    use_hq: bool,
    top_pvalues: usize,
) -> Option<Vec<f64>> {
    let scores = wait_for(id, || SCORES.get()).await;

    if top_pvalues == 0 || edit_count == 0 {
        let mut data = Vec::new();
        for (i, s) in scores.scores().enumerate() {
            data.push(if use_hq {
                s.edited_score_hq(edit_count)?
            } else {
                s.edited_score(edit_count)?
            });

            if i % 400 == 0 {
                yield_now(id, Some(i)).await;
            }
        }
        Some(data)
    } else {
        let mut top_variants = top_pvalue_variants(id, use_hq).await;
        if top_variants.len() < top_pvalues {
            return None;
        }
        top_variants = &top_variants[0..top_pvalues];
        if top_variants.len() < edit_count.unsigned_abs() {
            return None;
        }

        yield_now(id, None).await;

        let edits = scores.scores.keys().map(|id| {
            get_top_edits(top_variants, id, use_hq, edit_count)
                .into_iter()
                .map(|(s, g)| actual_edit(s, g, use_hq, edit_count).unwrap())
                .sum::<NotNan<f64>>()
        });

        let base_scores = Box::pin(compute_scores_(id, 0, use_hq, 0)).await?;

        Some(
            futures::stream::iter(base_scores.into_iter().zip(edits))
                .enumerate()
                .then(|(i, (base, edit))| async move {
                    if i % 400 == 0 {
                        yield_now(id, Some(i)).await;
                    }
                    base + *edit
                })
                .collect()
                .await,
        )
    }
}

/// `async` to give it a chance to be interrupted.
async fn compute_edit_analysis(
    id: HandlerId,
    edit_count: isize,
    use_hq: bool,
    top_pvalues: usize,
) -> Option<OutputGetEditAnalysis> {
    log::info!("[Worker][{id:?}] Computing edit analysis.");

    let scores = wait_for(id, || SCORES.get()).await;

    if top_pvalues == 0 || edit_count == 0 {
        let decreasing = edit_count < 0;
        let edit_count = edit_count.unsigned_abs();

        let mut edits = BTreeMap::new();
        for (i, score) in scores.scores().enumerate() {
            let all_edits = match (use_hq, decreasing) {
                (true, true) => &score.worst_hq,
                (true, false) => &score.best_hq,
                (false, true) => &score.worst,
                (false, false) => &score.best,
            };
            let all_edits = all_edits.get(..edit_count)?;
            for variant in all_edits.iter() {
                let key = variant.key.clone().key;
                let edit = edits
                    .entry(key.clone())
                    .or_insert_with(|| VariantInfo::new(key, scores));

                edit.count += 1;
                edit.ploidy_dosage.increment(VariantSampleInfo {
                    dosage: variant.dosage,
                    ploidy: variant.ploidy,
                });
            }

            if i % 400 == 0 {
                yield_now(id, Some(i)).await;
            }
        }

        Some(OutputGetEditAnalysis {
            use_hq,
            edits: edits.into_values().collect(),
        })
    } else {
        let mut top_variants = top_pvalue_variants(id, use_hq).await;
        if top_variants.len() < top_pvalues {
            return None;
        }
        top_variants = &top_variants[0..top_pvalues];
        if top_variants.len() < edit_count.unsigned_abs() {
            return None;
        }

        yield_now(id, None).await;

        let mut edits = BTreeMap::new();
        for (i, sample) in scores.scores.keys().enumerate() {
            for (s, g) in get_top_edits(top_variants, sample, use_hq, edit_count) {
                let key = SummaryStatKey::new(s);
                let edit = edits.entry(key).or_insert_with(|| VariantInfo {
                    count: 0,
                    stat: AsJson(s.clone()),
                    ploidy_dosage: CountingBTreeSet::new(),
                });
                edit.count += 1;
                edit.ploidy_dosage.increment(VariantSampleInfo {
                    dosage: g.dosage,
                    ploidy: g.ploidy,
                });
            }
            if i % 400 == 0 {
                yield_now(id, Some(i)).await;
            }
        }

        Some(OutputGetEditAnalysis {
            use_hq,
            edits: edits.into_values().collect(),
        })
    }
}

/// The variants with the largest effects in the slice.
fn get_top_edits<'a>(
    top_variants: &'a [&TopPValueVariant],
    id: &str,
    use_hq: bool,
    edit_count: isize,
) -> Vec<(&'a SummaryStats<GRCh38Contig>, SampleGenotype)> {
    let positive_edits = edit_count >= 0;

    let mut edits: Vec<_> = top_variants
        .iter()
        .map(|v| (&v.stat, v.genotypes[id]))
        .collect();

    edits.sort_unstable_by_key(|(s, g)| {
        let edit = actual_edit(s, *g, use_hq, edit_count).unwrap();
        if positive_edits {
            -edit
        } else {
            edit
        }
    });

    edits.truncate(edit_count.unsigned_abs());

    edits
}

fn actual_edit(
    s: &SummaryStats<GRCh38Contig>,
    g: SampleGenotype,
    use_hq: bool,
    edit_count: isize,
) -> Option<NotNan<f64>> {
    match (use_hq, edit_count >= 0) {
        (true, true) => s.max_edit_hq(g.dosage, g.ploidy),
        (true, false) => s.min_edit_hq(g.dosage, g.ploidy),
        (false, true) => s.max_edit(g.dosage, g.ploidy),
        (false, false) => s.min_edit(g.dosage, g.ploidy),
    }
}

async fn wait_for<T>(id: HandlerId, f: impl Fn() -> Option<T>) -> T {
    loop {
        if let Some(t) = f() {
            return t;
        }

        log::info!("[Worker][{id:?}] Data not ready.");
        leptos_ext::util::sleep(Duration::from_millis(200)).await;
    }
}

/// Gives a chance for new requests to replace the current one.
async fn yield_now(id: HandlerId, i: Option<usize>) {
    log::info!("[Worker][{id:?}] Yielding: {i:?}");
    utile::time::sleep(Duration::from_millis(0)).await;
}

async fn fetch_scores(id: HandlerId, origin: &str, file: String) -> io::Result<Scores> {
    let url = format!("{origin}/data/pan_ukbb/scores/{file}.json.br");

    log::info!("[Worker][{id:?}] Loading initial data from URL: {url}");

    let mut scores: Scores = load_large_json(&url).await?;

    // "Based on these thresholds, two individuals, NA21310 and HG02300,
    // were listed as males, but had genotypes consistent with females"
    // https://www.biorxiv.org/content/10.1101/078600v1.full.pdf
    scores.scores.remove("HG02300");
    scores.scores.remove("NA21310");

    log::info!("[Worker][{id:?}] Data loaded from {url}");

    Ok(scores)
}
async fn fetch_pvalues(
    id: HandlerId,
    origin: &str,
    file: String,
) -> io::Result<PhenotypeTopPValues> {
    let url = format!("{origin}/data/pan_ukbb/pvalues/{file}.json.br");

    log::info!("[Worker][{id:?}] Loading initial data from URL: {url}");

    let mut pvalues: PhenotypeTopPValues = load_large_json(&url).await?;

    // "Based on these thresholds, two individuals, NA21310 and HG02300,
    // were listed as males, but had genotypes consistent with females"
    // https://www.biorxiv.org/content/10.1101/078600v1.full.pdf
    pvalues.samples.remove("HG02300");
    pvalues.samples.remove("NA21310");

    for variant in pvalues.top_variants.values_mut() {
        variant.genotypes.remove("HG02300");
        variant.genotypes.remove("NA21310");
    }

    log::info!("[Worker][{id:?}] Data loaded from {url}");

    Ok(pvalues)
}
pub async fn fetch_manifest(origin: String) -> io::Result<Vec<PhenotypeManifestEntry>> {
    let url = format!("{origin}/data/pan_ukbb/phenotype_manifest.tsv");

    log::info!("[App] Fetching manifest from: {url}");

    let entries = PhenotypeManifestEntry::load_async(UrlResource::new(&url)?).await?;

    log::info!("[App] Manifest loaded.");

    Ok(entries)
}

/// Hack that loads the (compressed!) raw data asyncronously, then deserializes it in-place.
/// Avoids very large allocation of uncompressed data.
// TODO: move resource to uitle.
async fn load_large_json<T: DeserializeOwned>(url: &str) -> io::Result<T> {
    struct InMemoryResource {
        data: RefCell<Option<Vec<u8>>>,
    }
    impl RawResource for InMemoryResource {
        const NAMESPACE: &'static str = "memory";
        fn key(&self) -> String {
            unreachable!()
        }
        fn compression(&self) -> Option<Compression> {
            None
        }

        type Reader = std::io::Cursor<Vec<u8>>;
        fn size(&self) -> std::io::Result<u64> {
            unreachable!()
        }
        fn read(&self) -> std::io::Result<Self::Reader> {
            Ok(std::io::Cursor::new(self.data.borrow_mut().take().unwrap()))
        }

        type AsyncReader = std::io::Cursor<Vec<u8>>;
        async fn size_async(&self) -> std::io::Result<u64> {
            unreachable!()
        }
        async fn read_async(&self) -> std::io::Result<Self::AsyncReader> {
            unreachable!()
        }
    }

    let raw = UrlResource::new(url)?.read_vec_async().await?;

    InMemoryResource {
        data: RefCell::new(Some(raw)),
    }
    .decompressed_with(Compression::Brotli)
    .read_json()
}

fn spawn_task(
    task: &'static LocalKey<RefCell<Option<(HandlerIdOrd, AbortHandle)>>>,
    id: HandlerId,
    f: impl Future<Output = ()> + 'static,
) {
    let id = HandlerIdOrd(id);
    let (abort_handle, abort_registration) = AbortHandle::new_pair();
    let f = Abortable::new(f, abort_registration);
    let f = async move {
        let _ = f.await;
    };

    task.with(move |task| {
        let task = &mut *task.borrow_mut();
        *task = Some(match task {
            Some((current_id, current_abort_handle)) if *current_id < id => {
                log::info!("[Worker][{id:?}] Aborting {current_id:?}.");
                current_abort_handle.abort();
                wasm_bindgen_futures::spawn_local(f);
                (id, abort_handle)
            }
            Some(_) => {
                log::info!("[Worker][{id:?}] Stale request, dropping.");
                return;
            }
            None => {
                log::info!("[Worker][{id:?}] First request.");
                wasm_bindgen_futures::spawn_local(f);
                (id, abort_handle)
            }
        });
    });
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct HandlerIdOrd(HandlerId);
impl std::fmt::Debug for HandlerIdOrd {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        <HandlerId as std::fmt::Debug>::fmt(&self.0, f)
    }
}
impl HandlerIdOrd {
    pub fn raw(self) -> usize {
        format!("{self:?}")
            .strip_prefix("HandlerId(")
            .unwrap()
            .strip_suffix(")")
            .unwrap()
            .parse()
            .unwrap()
    }
}
impl PartialOrd for HandlerIdOrd {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HandlerIdOrd {
    fn cmp(&self, other: &Self) -> Ordering {
        self.raw().cmp(&other.raw())
    }
}

// Horrible hack
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct AsJson<T>(T);
impl<T> AsJson<T> {
    pub fn inner(self) -> T {
        self.0
    }
}
impl<T> Deref for AsJson<T> {
    type Target = T;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl<T> DerefMut for AsJson<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
impl<T> Serialize for AsJson<T>
where
    T: Serialize,
{
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serde_json::to_string(&self.0)
            .map_err(serde::ser::Error::custom)?
            .serialize(serializer)
    }
}
impl<'de, T> Deserialize<'de> for AsJson<T>
where
    T: Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s: &str = serde::Deserialize::deserialize(deserializer)?;
        let v = serde_json::from_str(s).map_err(serde::de::Error::custom)?;
        Ok(Self(v))
    }
}
