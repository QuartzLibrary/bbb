#![feature(btree_set_entry)]

pub mod send_option;

use futures::stream::{AbortHandle, Abortable};
use gloo_worker::{HandlerId, Worker, WorkerScope};
use ordered_float::NotNan;
use serde::{Deserialize, Serialize};
use std::{
    cell::RefCell,
    cmp::Ordering,
    future::Future,
    io,
    ops::{Deref, DerefMut},
    rc::Rc,
    time::Duration,
};

use pan_ukbb::PhenotypeManifestEntry;
use utile::resource::{Compression, RawResourceExt, UrlResource};

use analysis::scores::{Scores, Stats};

thread_local! {
    static GLOBAL_STATE: RefCell<GlobalState> = const { RefCell::new(GlobalState::None) };
    static TASK: RefCell<Option<(HandlerIdOrd, AbortHandle)>> = const { RefCell::new(None) };
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum GlobalState {
    None,
    Loading {
        file: String,
    },
    Loaded {
        file: String,
        scores: Rc<Scores>,
        normalised: Rc<Scores>,
    },
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
    },
}

#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
#[expect(clippy::large_enum_variant)]
pub enum Output {
    Init(OutputInit),
    GetScores(Option<OutputGetScores>),
}
#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
pub struct OutputInit {
    // Horrible hack (AsJson)
    pub phenotype: AsJson<PhenotypeManifestEntry>,

    pub stats: Option<Stats>,
    pub stats_hq: Option<Stats>,

    pub norm_stats: Option<Stats>,
    pub norm_stats_hq: Option<Stats>,
}
#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
pub struct OutputGetScores {
    pub scores: Vec<f64>,

    pub stats: Stats,
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
                GLOBAL_STATE.with(|f| {
                    assert_eq!(&*f.borrow(), &GlobalState::None);
                    f.replace(GlobalState::Loading { file: file.clone() });
                });
                let scope = scope.clone();
                wasm_bindgen_futures::spawn_local(async move {
                    let scores = fetch_scores(id, origin, file.clone()).await.unwrap();
                    let normalised = scores.normalised();
                    scope.respond(
                        id,
                        Output::Init(OutputInit {
                            phenotype: AsJson(scores.phenotype.clone()),

                            stats: scores.stats(),
                            stats_hq: scores.stats_hq(),

                            norm_stats: normalised.stats(),
                            norm_stats_hq: normalised.stats_hq(),
                        }),
                    );

                    GLOBAL_STATE.with(|f| {
                        f.replace(GlobalState::Loaded {
                            file,
                            scores: Rc::new(scores),
                            normalised: Rc::new(normalised),
                        });
                    });
                });
            }
            Input::GetScores {
                edit_count,
                normalise,
                use_hq,
            } => {
                let scope = scope.clone();
                let f = async move {
                    loop {
                        if let Some(data) =
                            try_compute_scores(id, edit_count, normalise, use_hq).await
                        {
                            log::info!("[Worker][{id:?}] Sending scores.");
                            scope.respond(id, data);
                            break;
                        }

                        log::info!("[Worker][{id:?}] Data not ready.");

                        leptos_ext::util::sleep(Duration::from_millis(200)).await;
                    }
                };
                spawn_task(id, f);
            }
        }
    }
}

fn spawn_task(id: HandlerId, f: impl Future<Output = ()> + 'static) {
    let id = HandlerIdOrd(id);
    let (abort_handle, abort_registration) = AbortHandle::new_pair();
    let f = Abortable::new(f, abort_registration);
    let f = async move {
        let _ = f.await;
    };

    TASK.with(move |task| {
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

/// `async` to give it a chance to be interrupted.
async fn try_compute_scores(
    id: HandlerId,
    edit_count: isize,
    normalise: bool,
    use_hq: bool,
) -> Option<Output> {
    let scores = GLOBAL_STATE.with(|f| {
        let state = &mut *f.borrow_mut();
        match state {
            GlobalState::Loaded {
                file: _,
                scores: _,
                normalised,
            } if normalise => Some(normalised.clone()),
            GlobalState::Loaded {
                file: _,
                scores,
                normalised: _,
            } => Some(scores.clone()),
            GlobalState::None | GlobalState::Loading { .. } => None,
        }
    })?;

    log::info!("[Worker][{id:?}] Computing scores.");

    let Some(values) = compute_scores(edit_count, use_hq, &scores).await else {
        return Some(Output::GetScores(None));
    };

    let iter = || values.iter().copied();
    let mean = analysis::util::mean(iter()).unwrap();
    let std_dev = analysis::util::std_dev(mean, iter()).unwrap();
    let iter = || values.iter().copied().map(|v| NotNan::new(v).unwrap());
    let min = *iter().min().unwrap();
    let max = *iter().max().unwrap();

    Some(Output::GetScores(Some(OutputGetScores {
        scores: values,

        stats: Stats {
            mean,
            std_dev,
            min,
            max,
        },
    })))
}
/// `async` to give it a chance to be interrupted.
async fn compute_scores(edit_count: isize, use_hq: bool, scores: &Scores) -> Option<Vec<f64>> {
    let mut data = Vec::new();
    for (i, s) in scores.scores().enumerate() {
        data.push(if use_hq {
            s.edited_score_hq(edit_count)?
        } else {
            s.edited_score(edit_count)?
        });

        if i & 100 == 0 {
            yield_now().await;
        }
    }
    Some(data)
}

async fn yield_now() {
    utile::time::sleep(Duration::from_millis(0)).await;
}

async fn fetch_scores(id: HandlerId, origin: String, file: String) -> io::Result<Scores> {
    let url = format!("{origin}/{file}");

    log::info!("[Worker][{id:?}] Loading initial data from URL: {url}");

    let mut scores: Scores = UrlResource::new(&url)?
        .decompressed_with(Compression::Gzip)
        .read_json_async()
        .await?;

    // TODO: hack, why is this sample so weird?
    scores.scores.remove("HG02300");

    log::info!("[Worker][{id:?}] Data loaded from {url}");

    Ok(scores)
}
pub async fn fetch_manifest(origin: String) -> io::Result<Vec<PhenotypeManifestEntry>> {
    let url = format!("{origin}/data/pan_ukbb/phenotype_manifest.tsv");

    log::info!("[App] Fetching manifest from: {url}");

    let entries = PhenotypeManifestEntry::load_async(UrlResource::new(&url)?).await?;

    log::info!("[App] Manifest loaded.");

    Ok(entries)
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
