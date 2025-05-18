use gloo_worker::Spawnable;
use leptos::{
    attr::Attribute,
    ev,
    html::{self},
    prelude::*,
    view, IntoView,
};
use leptos_ext::signal::{Load, ReadSignalExt, WriteSignalExt};
use plotly::{
    layout::{Shape, ShapeLine},
    Plot,
};
use std::{
    sync::{Arc, Mutex},
    time::Duration,
};
use thaw::*;

use pan_ukbb::{PhenotypeManifestEntry, Population};

use analysis::scores::Stats;

use bbb::{send_option::SendOption, Input, Output, OutputGetScores, OutputInit, WorkerStruct};

const BINS: usize = 100;

pub fn score_page(file: String) -> impl IntoView {
    let initial_state = ArcRwSignal::new(None);

    let bridge = WorkerStruct::spawner()
        .callback({
            let initial_state = initial_state.clone();
            move |output| match output {
                Output::Init(init) => initial_state.set(Some(init)),
                Output::GetScores(_) => unreachable!(),
            }
        })
        .spawn("/worker.js");

    bridge.send(Input::Set {
        file: format!("data/pan_ukbb/scores/{file}.json.bz"),
        origin: leptos::prelude::window().location().origin().unwrap(),
    });

    let edit_count: RwSignal<isize> = RwSignal::new(0);
    let use_hq_scores = RwSignal::new(false);
    let normalise = RwSignal::new(false);
    let show_stats = RwSignal::new(false);

    let bridge = Arc::new(SendOption::new_local(Some(bridge)));

    let data = RwSignal::new(Load::Loading);

    Signal::derive(move || Input::GetScores {
        edit_count: edit_count.get(),
        normalise: normalise.get(),
        use_hq: use_hq_scores.get(),
    })
    .rate_limit_leading(Duration::from_millis(10))
    .map_async(move |input| {
        let input = input.clone();
        let bridge = bridge.clone();
        async move {
            let (tx, rx) = futures::channel::oneshot::channel::<Output>();
            let bridge = {
                let callback = Mutex::new(Some(move |output| {
                    let _ = tx.send(output);
                }));
                SendOption::new_local(Some((*bridge).as_ref().unwrap().fork(Some(
                    move |output| {
                        log::info!("[App] Received output");
                        callback.lock().unwrap().take().unwrap()(output);
                    },
                ))))
            };
            log::info!("[App] Requesting {input:?}");
            bridge.as_ref().unwrap().send(input);
            let output = rx.await.unwrap();
            log::info!("[App] Received output (async)");
            match output {
                Output::Init(_init) => unreachable!(),
                Output::GetScores(scores) => scores,
            }
        }
    })
    .for_each_immediate({
        let mut score_signal: Option<ArcRwSignal<_>> = None;

        move |scores| match scores {
            Load::Loading => data.set(Load::Loading),
            Load::Ready(None) => data.set(Load::Ready(None)),
            Load::Ready(Some(scores)) => {
                let scores = scores.clone();
                let s = match score_signal.clone() {
                    Some(s) => {
                        s.set(scores);
                        s
                    }
                    None => ArcRwSignal::new(scores),
                };
                score_signal = Some(s.clone());
                if !data.read_untracked().is_ready() {
                    data.set(Load::Ready(Some(s.clone())));
                }
            }
        }
    });
    let loading = data.map(|o| !o.is_ready());
    let data = data.keep_if(|o| o.is_ready());

    move || {
        let initial_state = initial_state.get();

        let slider = view! { <Slider value={RwSignal::from(edit_count.double_bind(|i| *i as f64, |i| *i as isize))} min=-500. max=500.0 step=1.0 show_stops=false style="width: 100%;" />};
        let plus = html::button()
            .class("btn")
            .on(ev::click, move |_| {
                edit_count.update(|c| {
                    if *c < 500 - 1 {
                        *c += 1
                    }
                })
            })
            .child("+");
        let minus = html::button()
            .class("btn")
            .on(ev::click, move |_| {
                edit_count.update(|c| {
                    if *c > -500 + 1 {
                        *c -= 1
                    }
                })
            })
            .child("-");

        let data = data.read();
        let data = data.as_ref().cloned();

        html::div().class("score-viewer").child((
            html::h1().child(
                initial_state
                    .as_ref()
                    .map(|init| init.phenotype.description.clone())
                    .unwrap_or("Score Viewer".to_owned()),
            ),
            html::div()
                .class("chart-container")
                .class(("loading", loading))
                .child(match data.clone() {
                    Load::Ready(Some(data)) => render_chart(
                        data.into(),
                        initial_state.clone().unwrap(),
                        use_hq_scores.into(),
                        normalise.into(),
                        show_stats.into(),
                    )
                    .into_any(),
                    Load::Ready(None) => html::div()
                        .class("loading-text")
                        .child("Too many edits, no scores available.")
                        .into_any(),
                    Load::Loading => html::div()
                        .class("loading-text")
                        .child("Loading...")
                        .into_any(),
                }),
            html::div().class("controls").child((
                html::div().class("control-group").child((
                    html::span()
                        .class("control-label")
                        .child("Number of edits: "),
                    slider,
                    html::div()
                        .class("increment-buttons")
                        .child((minus, edit_count, plus)),
                )),
                html::div().class("control-group").child((
                    html::span().class("control-label").child("Use HQ scores: "),
                    view! { <Switch checked=use_hq_scores /> },
                )),
                html::div().class("control-group").child((
                    html::span().class("control-label").child("Normalise: "),
                    view! { <Switch checked=normalise /> },
                )),
                html::div().class("control-group").child((
                    html::span().class("control-label").child("Show stats: "),
                    view! { <Switch checked=show_stats /> },
                )),
            )),
            move || {
                initial_state.as_ref().map(|init| {
                    render_full_phenotype_info(
                        init,
                        data.clone(),
                        use_hq_scores.get(),
                        normalise.get(),
                    )
                })
            },
        ))
    }
}

fn render_full_phenotype_info(
    init: &OutputInit,
    data: Load<Option<ArcRwSignal<OutputGetScores>>>,
    use_hq: bool,
    normalise: bool,
) -> impl IntoView {
    fn label_value(label: &str, value: impl IntoView) -> impl IntoView {
        html::div().class("info-item").child((
            html::span().class("info-label").child(label.to_string()),
            html::span().class("info-value").child(value),
        ))
    }
    fn title(title: &'static str) -> impl IntoView {
        html::h2().class("section-title").child(title)
    }
    fn render_stats_section(t: &'static str, stats: Option<Stats>) -> impl IntoView {
        html::div().class("info-section").child((
            title(t),
            html::div().class("info-grid").child(
                stats
                    .map(|stats| {
                        (
                            html::div().class("stats-row").child((
                                label_value("Mean: ", format!("{:.4}", stats.mean)),
                                label_value("Std Dev: ", format!("{:.4}", stats.std_dev)),
                            )),
                            html::div().class("stats-row").child((
                                label_value("Min: ", format!("{:.4}", stats.min)),
                                label_value("Max: ", format!("{:.4}", stats.max)),
                            )),
                        )
                    })
                    .unwrap_or_else(|| {
                        (
                            html::div().class("stats-row").child((
                                label_value("Mean: ", "Loading...".to_owned()),
                                label_value("Std Dev: ", "Loading...".to_owned()),
                            )),
                            html::div().class("stats-row").child((
                                label_value("Min: ", "Loading...".to_owned()),
                                label_value("Max: ", "Loading...".to_owned()),
                            )),
                        )
                    })
                    .into_any(),
            ),
        ))
    }

    let PhenotypeManifestEntry {
        trait_type,
        phenocode,
        pheno_sex,
        coding,
        modifier,

        description,
        description_more,
        coding_description,
        category,
        in_max_independent_set,

        n_cases_full_cohort_both_sexes,
        n_cases_full_cohort_females,
        n_cases_full_cohort_males,
        n_cases_hq_cohort_both_sexes,
        n_cases_hq_cohort_females,
        n_cases_hq_cohort_males,

        pops,
        num_pops: _,
        pops_pass_qc,
        num_pops_pass_qc: _,

        n_cases_AFR,
        n_cases_AMR,
        n_cases_CSA,
        n_cases_EAS,
        n_cases_EUR,
        n_cases_MID,

        n_controls_AFR,
        n_controls_AMR,
        n_controls_CSA,
        n_controls_EAS,
        n_controls_EUR,
        n_controls_MID,

        rhemc_25bin_50rv_h2_observed_AFR,
        rhemc_25bin_50rv_h2_observed_AMR,
        rhemc_25bin_50rv_h2_observed_CSA,
        rhemc_25bin_50rv_h2_observed_EAS,
        sldsc_25bin_h2_observed_EUR,
        rhemc_25bin_50rv_h2_observed_MID,

        rhemc_25bin_50rv_h2_observed_se_AFR,
        rhemc_25bin_50rv_h2_observed_se_AMR,
        rhemc_25bin_50rv_h2_observed_se_CSA,
        rhemc_25bin_50rv_h2_observed_se_EAS,
        sldsc_25bin_h2_observed_se_EUR,
        rhemc_25bin_50rv_h2_observed_se_MID,

        rhemc_25bin_50rv_h2_liability_AFR: _,
        rhemc_25bin_50rv_h2_liability_AMR: _,
        rhemc_25bin_50rv_h2_liability_CSA: _,
        rhemc_25bin_50rv_h2_liability_EAS: _,
        sldsc_25bin_h2_liability_EUR: _,
        rhemc_25bin_50rv_h2_liability_MID: _,

        rhemc_25bin_50rv_h2_liability_se_AFR: _,
        rhemc_25bin_50rv_h2_liability_se_AMR: _,
        rhemc_25bin_50rv_h2_liability_se_CSA: _,
        rhemc_25bin_50rv_h2_liability_se_EAS: _,
        sldsc_25bin_h2_liability_se_EUR: _,
        rhemc_25bin_50rv_h2_liability_se_MID: _,

        rhemc_25bin_50rv_h2_z_AFR: _,
        rhemc_25bin_50rv_h2_z_AMR: _,
        rhemc_25bin_50rv_h2_z_CSA: _,
        rhemc_25bin_50rv_h2_z_EAS: _,
        sldsc_25bin_h2_z_EUR: _,
        rhemc_25bin_50rv_h2_z_MID: _,

        lambda_gc_AFR: _,
        lambda_gc_AMR: _,
        lambda_gc_CSA: _,
        lambda_gc_EAS: _,
        lambda_gc_EUR: _,
        lambda_gc_MID: _,

        phenotype_qc_AFR: _,
        phenotype_qc_AMR: _,
        phenotype_qc_CSA: _,
        phenotype_qc_EAS: _,
        phenotype_qc_EUR: _,
        phenotype_qc_MID: _,

        filename: _,
        filename_tabix: _,
        aws_path: _,
        aws_path_tabix: _,
        md5_hex,
        size_in_bytes,
        md5_hex_tabix: _,
        size_in_bytes_tabix: _,
    } = init.phenotype.clone().inner();

    html::div().class("phenotype-info").child((
        // Basic Information
        html::div().class("info-section").child((
            title("Basic Information"),
            html::div().class("info-grid").child((
                label_value("Type: ", format!("{trait_type}")),
                label_value("Code: ", phenocode.clone()),
                label_value("Sex: ", format!("{pheno_sex}")),
                coding.map(|coding| label_value("Coding: ", coding)),
                modifier.map(|modifier| label_value("Modifier: ", format!("{modifier}"))),
                category.map(|category| label_value("Category: ", category)),
            )),
        )),
        // Description
        html::div().class("info-section").child((
            title("Description"),
            html::div().class("description-content").child((
                html::p()
                    .class("main-description")
                    .child(description.clone()),
                description_more.map(|d| html::p().class("detailed-description").child(d.clone())),
                coding_description.map(|d| {
                    html::div().class("coding-description").child((
                        html::h3()
                            .class("subsection-title")
                            .child("Coding Description"),
                        html::p().child(d.clone()),
                    ))
                }),
            )),
        )),
        // Statistics
        {
            let (title, stats) = match (use_hq, normalise) {
                (true, true) => ("Normalised Statistics (HQ)", init.norm_stats_hq),
                (true, false) => ("Original Statistics (HQ)", init.stats_hq),
                (false, true) => ("Normalised Statistics", init.norm_stats),
                (false, false) => ("Original Statistics", init.stats),
            };
            render_stats_section(title, stats)
        },
        move || match &data {
            Load::Ready(Some(data)) => {
                render_stats_section("Edited Statistics", Some(data.get().stats))
            }
            _ => render_stats_section("Edited Statistics", None),
        },
        // Case Counts
        html::div().class("info-section").child((
            title("Case Counts"),
            html::div().class("info-grid").child((
                label_value(
                    "Total Cases: ",
                    format!(
                        "{n_cases_full_cohort_both_sexes} (HQ: {})",
                        n_cases_hq_cohort_both_sexes.unwrap_or(0)
                    ),
                ),
                label_value(
                    "Female Cases: ",
                    format!(
                        "{n_cases_full_cohort_females} (HQ: {})",
                        n_cases_hq_cohort_females.unwrap_or(0)
                    ),
                ),
                label_value(
                    "Male Cases: ",
                    format!(
                        "{n_cases_full_cohort_males} (HQ: {})",
                        n_cases_hq_cohort_males.unwrap_or(0)
                    ),
                ),
            )),
        )),
        // Population Information
        html::div().class("info-section").child((
            title("Population Information"),
            html::div().class("info-grid").child((
                label_value(
                    "Populations: ",
                    format!(
                        "{} ({} pass QC)",
                        pops.iter()
                            .map(|p| p.to_string())
                            .collect::<Vec<_>>()
                            .join(", "),
                        pops_pass_qc
                            .iter()
                            .map(|p| p.to_string())
                            .collect::<Vec<_>>()
                            .join(", ")
                    ),
                ),
                label_value(
                    "In Max Independent Set: ",
                    if in_max_independent_set { "Yes" } else { "No" },
                ),
            )),
        )),
        // Population-specific Data
        html::div().class("info-section").child((
            title("Population-specific Data"),
            html::div().class("population-grid").child(
                Population::all()
                    .into_iter()
                    .map(|pop| {
                        let cases = match pop {
                            Population::Afr => n_cases_AFR,
                            Population::Amr => n_cases_AMR,
                            Population::Csa => n_cases_CSA,
                            Population::Eas => n_cases_EAS,
                            Population::Eur => n_cases_EUR,
                            Population::Mid => n_cases_MID,
                        };
                        let controls = match pop {
                            Population::Afr => n_controls_AFR,
                            Population::Amr => n_controls_AMR,
                            Population::Csa => n_controls_CSA,
                            Population::Eas => n_controls_EAS,
                            Population::Eur => n_controls_EUR,
                            Population::Mid => n_controls_MID,
                        };
                        let h2 = match pop {
                            Population::Afr => rhemc_25bin_50rv_h2_observed_AFR,
                            Population::Amr => rhemc_25bin_50rv_h2_observed_AMR,
                            Population::Csa => rhemc_25bin_50rv_h2_observed_CSA,
                            Population::Eas => rhemc_25bin_50rv_h2_observed_EAS,
                            Population::Eur => sldsc_25bin_h2_observed_EUR,
                            Population::Mid => rhemc_25bin_50rv_h2_observed_MID,
                        };
                        let h2_se = match pop {
                            Population::Afr => rhemc_25bin_50rv_h2_observed_se_AFR,
                            Population::Amr => rhemc_25bin_50rv_h2_observed_se_AMR,
                            Population::Csa => rhemc_25bin_50rv_h2_observed_se_CSA,
                            Population::Eas => rhemc_25bin_50rv_h2_observed_se_EAS,
                            Population::Eur => sldsc_25bin_h2_observed_se_EUR,
                            Population::Mid => rhemc_25bin_50rv_h2_observed_se_MID,
                        };

                        html::div().class("population-item").child((
                            html::h3().class("population-title").child(pop.to_string()),
                            html::div().class("population-details").child((
                                cases.map(|n| {
                                    html::div()
                                        .class("population-stat")
                                        .child(format!("Cases: {n}"))
                                }),
                                controls.map(|n| {
                                    html::div()
                                        .class("population-stat")
                                        .child(format!("Controls: {n}"))
                                }),
                                h2.map(|h| {
                                    html::div()
                                        .class("population-stat")
                                        .child(format!("h²: {h}"))
                                }),
                                h2_se.map(|se| {
                                    html::div()
                                        .class("population-stat")
                                        .child(format!("h² SE: {se}"))
                                }),
                            )),
                        ))
                    })
                    .collect_view(),
            ),
        )),
        // File Information
        html::div().class("info-section").child((
            title("File Information"),
            html::div().class("info-grid").child((
                label_value(
                    "File Size: ",
                    format!("{:.1} MB", size_in_bytes as f64 / 1_000_000.0),
                ),
                label_value("MD5: ", md5_hex.clone()),
            )),
        )),
    ))
}

fn render_chart(
    data: Signal<OutputGetScores>,
    init: OutputInit,
    use_hq: Signal<bool>,
    normalise: Signal<bool>,
    show_stats: Signal<bool>,
) -> impl IntoView {
    const DEEP_SKY_BLUE: &str = "#00BFFF";
    const HOT_PINK: &str = "#FF69B4";
    const LIME_GREEN: &str = "#32CD32";
    const GOLD: &str = "#FFD700";

    fn stat_line(x: f64, color: &'static str, name: &'static str) -> Shape {
        Shape::new()
            .x0(x)
            .x1(x)
            .y0(0.0)
            .y1(1.0)
            .y_ref("paper")
            .line(ShapeLine::new().color(color).width(2.0))
            .name(name)
    }
    fn unedited_stat_line(x: f64, color: &'static str, name: &'static str) -> Shape {
        Shape::new()
            .x0(x)
            .x1(x)
            .y0(0.0)
            .y1(1.0)
            .y_ref("paper")
            .line(ShapeLine::new().color(format!("{color}30")).width(1.0))
            .name(name)
    }

    let chart = Signal::derive_local(move || {
        let OutputGetScores {
            scores,
            stats:
                Stats {
                    mean,
                    std_dev,
                    min,
                    max,
                },
        } = data.get();

        let mut plot = plotly::Plot::new();
        plot.add_trace(plotly::Histogram::new(scores).n_bins_x(BINS).name("Scores"));

        let current_shapes = [
            stat_line(mean, DEEP_SKY_BLUE, "Mean"), // Deep Sky Blue
            stat_line(mean + std_dev, HOT_PINK, "Mean + StdDev"), // Hot Pink
            stat_line(mean - std_dev, HOT_PINK, "Mean - StdDev"), // Hot Pink
            stat_line(min, LIME_GREEN, "Min"),      // Lime Green
            stat_line(max, GOLD, "Max"),            // Gold
        ];

        let unedited_shapes = match (use_hq.get(), normalise.get()) {
            (true, true) => init.norm_stats_hq,
            (true, false) => init.stats_hq,
            (false, true) => init.norm_stats,
            (false, false) => init.stats,
        }
        .map(|stats| {
            let Stats {
                mean,
                std_dev,
                min,
                max,
            } = stats;
            [
                unedited_stat_line(mean, DEEP_SKY_BLUE, "Mean"), // Deep Sky Blue
                unedited_stat_line(mean + std_dev, HOT_PINK, "Mean + StdDev"), // Hot Pink
                unedited_stat_line(mean - std_dev, HOT_PINK, "Mean - StdDev"), // Hot Pink
                unedited_stat_line(min, LIME_GREEN, "Min"),      // Lime Green
                unedited_stat_line(max, GOLD, "Max"),            // Gold
            ]
        });

        let shapes = if show_stats.get() {
            current_shapes
                .into_iter()
                .chain(unedited_shapes.into_iter().flatten())
                .collect()
        } else {
            vec![]
        };

        plot.set_layout(
            plotly::Layout::new()
                .template(&*plotly::layout::themes::PLOTLY_DARK)
                .shapes(shapes),
        );

        plot
    });

    render_plotly(chart)
}

fn render_plotly(
    plot: Signal<Plot, LocalStorage>,
) -> html::HtmlElement<html::Div, impl Attribute + Clone, impl RenderHtml + Clone> {
    let id: u128 = rand::random();

    Effect::new(move |_| {
        plot.with(|plot| {
            let plot = plot.clone();
            let plot_update_task = async move {
                log::info!("[App] Updating plot");
                plotly::bindings::react(&id.to_string(), &plot).await;
                log::info!("[App] Plot updated");
            };
            wasm_bindgen_futures::spawn_local(plot_update_task);
        });
    });

    html::div().id(id).style("width: 100%; height: 100%;")
}
