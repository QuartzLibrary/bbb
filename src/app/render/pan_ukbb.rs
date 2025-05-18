use gloo_worker::{Spawnable, WorkerBridge};
use leptos::{
    attr::Attribute,
    ev,
    html::{self},
    prelude::*,
    view, IntoView,
};
use leptos_ext::{
    signal::{Load, ReadSignalExt, WriteSignalExt},
    util::Task,
};
use ordered_float::NotNan;
use plotly::{
    configuration::DisplayModeBar,
    layout::{Shape, ShapeLine},
    Plot,
};
use std::{
    sync::{Arc, Mutex},
    time::Duration,
};
use thaw::*;
use web_sys;

use pan_ukbb::{PhenotypeManifestEntry, Population, SummaryStats};

use analysis::scores::Stats;

use edit::{
    send_option::SendOption, util::PLOTLY_THEME, Input, Output, OutputGetEditAnalysis,
    OutputGetScores, OutputInit, VariantInfo, VariantSampleInfo, WorkerStruct,
};

const BINS: usize = 100;

#[derive(Debug, Clone)]
struct PageState {
    _file: String,

    bridge: Arc<SendOption<WorkerBridge<WorkerStruct>>>,

    init: ArcRwSignal<Option<OutputInit>>,

    /// The number of edits to apply to the scores.
    /// Negative means reducing the scores.
    edit_count: RwSignal<isize>,
    /// Use the HQ fields from the PanUKBB data.
    use_hq: RwSignal<bool>,
    /// Normalise the scores and edits.
    normalise: RwSignal<bool>,
    /// Restrict the edits to values in the top n p-values.
    /// If 0, no restriction is applied.
    top_pvalues: RwSignal<usize>,

    show_stats: RwSignal<bool>,
    full_range: RwSignal<bool>,
}

pub fn score_page(file: String) -> impl IntoView {
    let state = PageState::new(file);

    let scores = state.scores_signal();
    let edits = state.edits_signal();

    let loading = scores.map(|o| !o.is_ready());
    let scores = scores.keep_if(|o| o.is_ready());

    move || {
        let PageState {
            _file: _,
            bridge: _,
            init,
            edit_count,
            use_hq,
            normalise,
            top_pvalues: _,
            show_stats: _,
            full_range: _,
        } = state.clone();

        let init = init.get();

        html::div().class("page-content score-viewer").child((
            html::a().href("/").class("back-link").child("< Index"),
            html::h1().child(
                init.as_ref()
                    .map(|init| init.phenotype.description.clone())
                    .unwrap_or("Loading…".to_owned()),
            ),
            init.as_ref()
                .and_then(|init| init.phenotype.description_more.clone())
                .map(|d| html::p().child(d)),
            html::div()
                .class("card card-chart")
                .class(("loading", loading))
                .child({
                    let state = state.clone();
                    move || match scores.get() {
                        Load::Ready(Some(scores)) => {
                            render_main_chart(&state, scores.into()).into_any()
                        }
                        Load::Ready(None) => html::div()
                            .class("loading-text")
                            .child("Too many edits, no scores available.")
                            .into_any(),
                        Load::Loading => html::div()
                            .class("loading-text")
                            .child("Loading…")
                            .into_any(),
                    }
                }),
            render_controls(&state),
            init.as_ref().map(|init| {
                let init = init.clone();
                move || {
                    render_full_phenotype_info(
                        &init,
                        scores.get(),
                        loading,
                        use_hq.get(),
                        normalise.get(),
                    )
                }
            }),
            init.as_ref()
                .map(|init| render_edit_analysis_panel(init.clone(), edits, edit_count.into())),
            html::p()
                .class("footnote")
                .child("All genomic positions are on GRCh38/hg38 and 1-based."),
        ))
    }
}

fn render_main_chart(state: &PageState, data: Signal<OutputGetScores>) -> impl IntoView {
    const DEEP_SKY_BLUE: &str = "#00BFFF";
    const HOT_PINK: &str = "#FF69B4";
    const LIME_GREEN: &str = "#32CD32";
    const GOLD: &str = "#FFD700";

    const MEAN: &str = "mean";
    const MEAN_PLUS_STDDEV: &str = "+sd";
    const MEAN_MINUS_STDDEV: &str = "-sd";
    const MIN: &str = "min";
    const MAX: &str = "max";

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

    fn create_stat_annotation(
        x: f64,
        text: &str,
        color: &'static str,
        opacity: f64,
    ) -> plotly::layout::Annotation {
        plotly::layout::Annotation::new()
            .x(x)
            .y(1.1)
            .y_ref("paper")
            .text(text)
            .show_arrow(false)
            .font(plotly::common::Font::new().color(color))
            .opacity(opacity)
    }

    let PageState {
        init,
        edit_count,
        use_hq,
        normalise,
        show_stats,
        full_range,
        ..
    } = state.clone();

    let init = init
        .get()
        .expect("initial state should be ready before scores");

    let zero_edits = edit_count.is(0).dedup();

    let chart = Signal::derive_local(move || {
        // Untracked to avoid re-rendering before the data is ready. (Always trigger a change in the data.)
        let use_hq = use_hq.get_untracked();
        let normalise = normalise.get_untracked();
        let zero_edits = zero_edits.get_untracked();

        let OutputGetScores { scores, stats } = data.get();
        let Stats {
            mean,
            std_dev,
            min,
            max,
        } = stats;

        let mut plot = plotly::Plot::new();
        plot.add_trace(plotly::Histogram::new(scores).n_bins_x(BINS).name("Scores"));

        let current_shapes = [
            stat_line(mean, DEEP_SKY_BLUE, MEAN),
            stat_line(mean + std_dev, HOT_PINK, MEAN_PLUS_STDDEV),
            stat_line(mean - std_dev, HOT_PINK, MEAN_MINUS_STDDEV),
            stat_line(min, LIME_GREEN, MIN),
            stat_line(max, GOLD, MAX),
        ];

        let unedited_shapes = match (use_hq, normalise) {
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
                unedited_stat_line(mean, DEEP_SKY_BLUE, MEAN),
                unedited_stat_line(mean + std_dev, HOT_PINK, MEAN_PLUS_STDDEV),
                unedited_stat_line(mean - std_dev, HOT_PINK, MEAN_MINUS_STDDEV),
                unedited_stat_line(min, LIME_GREEN, MIN),
                unedited_stat_line(max, GOLD, MAX),
            ]
        });

        let shapes = if show_stats.get() {
            let unedited_shapes = match unedited_shapes {
                Some(_) if zero_edits => None,
                Some(shapes) => Some(shapes),
                None => None,
            };

            current_shapes
                .into_iter()
                .chain(unedited_shapes.into_iter().flatten())
                .collect()
        } else {
            vec![]
        };

        // Create annotations for statistical lines
        let annotations = if show_stats.get() {
            vec![
                create_stat_annotation(mean, MEAN, DEEP_SKY_BLUE, 1.0),
                create_stat_annotation(mean + std_dev, MEAN_PLUS_STDDEV, HOT_PINK, 1.0),
                create_stat_annotation(mean - std_dev, MEAN_MINUS_STDDEV, HOT_PINK, 1.0),
                create_stat_annotation(min, MIN, LIME_GREEN, 1.0),
                create_stat_annotation(max, MAX, GOLD, 1.0),
            ]
        } else {
            vec![]
        };

        let mut layout = plotly::Layout::new()
            .template(PLOTLY_THEME.get())
            .shapes(shapes)
            .annotations(annotations)
            .margin(
                plotly::layout::Margin::new()
                    .top(60)
                    .right(60)
                    .bottom(60)
                    .left(60),
            );

        if full_range.get() {
            let range = match (use_hq, normalise) {
                (true, true) => init.norm_full_score_range_hq.clone(),
                (true, false) => init.full_score_range_hq.clone(),
                (false, true) => init.norm_full_score_range.clone(),
                (false, false) => init.full_score_range.clone(),
            };

            if let Some(range) = range {
                layout = layout.x_axis(
                    plotly::layout::Axis::new()
                        .range(vec![range.start, range.end])
                        .auto_range(false),
                );
            }
        }

        plot.set_layout(layout);

        plot
    });

    render_plotly(chart)
}

fn render_controls(state: &PageState) -> impl IntoView {
    const EDIT_COUNT_TIP: &str = "The number of edits to apply to the samples.";
    const TOP_PVALUES_TIP: &str = "Pick the edits only from the top-n variants with the lowest p-values. Set to 0 to pick from all variants.";

    let PageState {
        edit_count,
        use_hq,
        normalise,
        show_stats,
        full_range,
        top_pvalues,
        ..
    } = state.clone();

    let buttons = [
        (use_hq, "Use HQ", "Use the HQ fields from the PanUKBB data."),
        (normalise, "Normalise", "Normalise the scores to the range of the original data."),
        (show_stats, "Show stats", "Show some summary statistics in the chart above."),
        (full_range, "Full range", "Force the chart to show the full range of edits (makes it easert to eyeball the shift in the distribution)."),
    ];

    html::div().class("card controls").child((
        slider_with_controls(
            "Edits: ",
            edit_count
                .double_bind(
                    |i| *i as f64,
                    move |edit| {
                        let edit = *edit as isize;

                        // To make a certain number of edits, we need to pick from at least that many p-values.
                        let top = top_pvalues.get_untracked();
                        if top != 0 && top < edit.unsigned_abs() {
                            top_pvalues.set(edit.unsigned_abs());
                        }

                        edit
                    },
                )
                .into(),
            -500.0,
            500.0,
            1.0,
            "0",
            Some(EDIT_COUNT_TIP),
        ),
        slider_with_controls(
            "Top p-values: ",
            top_pvalues
                .double_bind(
                    |i| *i as f64,
                    move |top| {
                        let top = *top as usize;

                        let edit = edit_count.get_untracked();
                        if top != 0 && top < edit.unsigned_abs() {
                            edit_count.set(top as isize);
                        }

                        top
                    },
                )
                .into(),
            0.0,
            3_000.0,
            1.0,
            "any",
            Some(TOP_PVALUES_TIP),
        ),
        html::div()
            .class("controls-row")
            .child(buttons.map(|(signal, label, tip)| toggle_button(signal, label, Some(tip)))),
    ))
}

fn slider_with_controls(
    name: &'static str,
    signal: RwSignal<f64>,
    min: f64,
    max: f64,
    step: f64,
    zero: &'static str,
    tip: Option<&'static str>,
) -> impl IntoView {
    let slider = view! { <Slider value={signal} min={min} max={max} step={step} show_stops=false style="width: 100%;" />};
    let plus = html::button()
        .class("btn")
        .attr("disabled", signal.map(move |v| v == &max))
        .on(ev::click, move |_| {
            signal.update(|c| {
                if *c < max - step {
                    *c += step
                }
            })
        })
        .child("+");
    let minus = html::button()
        .class("btn")
        .attr("disabled", signal.map(move |v| v == &min))
        .on(ev::click, move |_| {
            signal.update(|c| {
                if *c > min + step {
                    *c -= step
                }
            })
        })
        .child("-");

    html::div()
        .class("control-edits")
        .class(("tooltip", tip.is_some()))
        .attr("data-text", tip)
        .child((
            html::span().class("control-label").child(name),
            slider,
            html::div().class("increment-buttons").child((
                minus,
                html::span()
                    .class("slider-value")
                    .child(signal.map(move |v| {
                        if *v == 0.0 {
                            zero.to_owned()
                        } else {
                            v.to_string()
                        }
                    })),
                plus,
            )),
        ))
}

fn render_full_phenotype_info(
    init: &OutputInit,
    data: Load<Option<ArcRwSignal<OutputGetScores>>>,
    loading: Signal<bool>,
    use_hq: bool,
    normalise: bool,
) -> impl IntoView {
    fn label_value(label: &str, value: impl IntoView) -> impl IntoView {
        html::div().class("info-item").child((
            html::span().class("info-label").child(label.to_string()),
            html::span().class("info-value").child(value),
        ))
    }
    fn render_stats_section(t: &'static str, stats: Stats, loading: Signal<bool>) -> impl IntoView {
        html::div()
            .class("card")
            .class(("loading", loading))
            .child((
                html::h2().child(t),
                html::div().class("info-grid").child((
                    html::div().class("stats-row").child((
                        label_value("Mean: ", format!("{:.4}", stats.mean)),
                        label_value("Std Dev: ", format!("{:.4}", stats.std_dev)),
                    )),
                    html::div().class("stats-row").child((
                        label_value("Min: ", format!("{:.4}", stats.min)),
                        label_value("Max: ", format!("{:.4}", stats.max)),
                    )),
                )),
            ))
    }
    fn render_stats_section_placeholder(
        t: &'static str,
        msg: &'static str,
        loading: Signal<bool>,
    ) -> impl IntoView {
        html::div()
            .class("card")
            .class(("loading", loading))
            .child((
                html::h2().child(t),
                html::div().class("info-grid").child((
                    html::div()
                        .class("stats-row")
                        .child((label_value("Mean: ", msg), label_value("Std Dev: ", msg))),
                    html::div()
                        .class("stats-row")
                        .child((label_value("Min: ", msg), label_value("Max: ", msg))),
                )),
            ))
    }

    let PhenotypeManifestEntry {
        trait_type,
        phenocode,
        pheno_sex,
        coding,
        modifier,

        description: _,
        description_more: _,
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
        md5_hex: _,
        size_in_bytes: _,
        md5_hex_tabix: _,
        size_in_bytes_tabix: _,
    } = init.phenotype.clone().inner();

    let phenocode_url = format!("https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id={phenocode}");

    let never_loading = Signal::derive(move || false);

    (
        // Statistics
        {
            let (title, stats) = match (use_hq, normalise) {
                (true, true) => ("Normalised Statistics (HQ)", init.norm_stats_hq),
                (true, false) => ("Original Statistics (HQ)", init.stats_hq),
                (false, true) => ("Normalised Statistics", init.norm_stats),
                (false, false) => ("Original Statistics", init.stats),
            };
            match stats {
                Some(stats) => render_stats_section(title, stats, never_loading).into_any(),
                None => render_stats_section_placeholder(title, "N/A", never_loading).into_any(),
            }
        },
        move || match &data {
            Load::Ready(Some(data)) => {
                render_stats_section("Edited Statistics", data.get().stats, loading).into_any()
            }
            Load::Ready(None) => {
                render_stats_section_placeholder("Edited Statistics", "N/A", loading).into_any()
            }
            Load::Loading => {
                render_stats_section_placeholder("Edited Statistics", "Loading…", loading)
                    .into_any()
            }
        },
        // Basic Information
        html::div().class("card").child((
            html::h2().child("Basic Information"),
            html::div().class("info-grid").child((
                label_value("Type: ", format!("{trait_type}")),
                label_value(
                    "Code: ",
                    html::a().href(phenocode_url).child(phenocode.clone()),
                ),
                label_value("Sex: ", format!("{pheno_sex}")),
                coding.map(|coding| label_value("Coding: ", coding)),
                coding_description.map(|coding_description| {
                    label_value("Coding Description: ", coding_description)
                }),
                modifier.map(|modifier| label_value("Modifier: ", format!("{modifier}"))),
                category.map(|category| label_value("Category: ", category)),
                label_value(
                    "Covariates: ",
                    "Age, Sex, Age * Sex, Age², Age² * Sex, and the first 10 PCs",
                ),
            )),
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
        // Population-specific Data
        html::div().class("card").child((
            html::h2().child("Population breakdown"),
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
                            html::h3().child(pop.to_string()),
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
    )
}

fn render_edit_analysis_panel(
    init: OutputInit,
    edits: Signal<Load<Option<ArcRwSignal<OutputGetEditAnalysis>>>>,
    edit_count: Signal<isize>,
) -> impl IntoView {
    let loading = edits.map(|o| !o.is_ready());
    let edits = edits.keep_if(|o| o.is_ready());

    let sort_by_pvalue = RwSignal::new(true);
    let sort_by_effect = sort_by_pvalue.double_bind(|v| !*v, |v| !*v);

    let sorted_edits = Signal::derive(move || {
        let mut data = edits.get().ready()??.get();
        let OutputGetEditAnalysis { use_hq, edits } = &mut data;
        if !sort_by_pvalue.get() {
            // They are pre-sorted by p-value.
            edits.sort_unstable_by_key(|a| {
                let beta = if *use_hq {
                    a.stat.beta_meta_hq.unwrap()
                } else {
                    a.stat.beta_meta.unwrap()
                };
                NotNan::new(-beta.abs()).unwrap()
            });
        }
        Some(data)
    });

    (
        html::div()
            .class("card card-chart")
            .class(("loading", loading))
            .child(move || match edits.get() {
                Load::Loading => html::div()
                    .class("loading-text")
                    .child("Loading…")
                    .into_any(),
                Load::Ready(None) => html::div()
                    .class("loading-text")
                    .child("Too many edits.")
                    .into_any(),
                Load::Ready(Some(_)) if edit_count.is(0).get() => html::div()
                    .class("loading-text")
                    .child("Apply at least one edit.")
                    .into_any(),
                Load::Ready(Some(analysis)) => {
                    render_edit_analysis_chart(edit_count, init.sample_count, analysis.into())
                        .into_any()
                }
            }),
        html::div().class("card sort-settings").child((
            "Sort by: ",
            [
                (sort_by_pvalue, "p-value"),
                (sort_by_effect.into(), "effect"),
            ]
            .map(|(signal, label)| toggle_button(signal, label, None)),
        )),
        move || {
            // Rendering long lists can be kind of slow, this is a horrible hack that just renders it in chunks,
            // breaking up the work and keeping the interface responsive rather than implementing a proper virtual list.
            // But also, why is it slow at all? 500 elements is not that much.

            let len = sorted_edits.with(|s| s.as_ref().map(|s| s.edits.len()).unwrap_or(0));

            let loading = {
                let l = ArcRwSignal::new(loading.get_untracked());
                loading.for_each_after_first({
                    let l = l.clone();
                    move |v| l.set(*v)
                });
                l
            };

            let pieces = (len / 20) + 1;

            (0..pieces)
                .map(|i| {
                    let start = (i * 20).min(len);
                    let end = (start + 20).min(len);

                    let loading = loading.clone();

                    lazy(
                        move || {
                            sorted_edits.with(|sorted_edits| {
                                let sorted_edits = sorted_edits.as_ref()?;
                                Some(
                                    sorted_edits.edits[start..end]
                                        .iter()
                                        .map(|variant_info| {
                                            edit_table_line(
                                                variant_info,
                                                sorted_edits.use_hq,
                                                loading.clone().into(),
                                            )
                                        })
                                        .collect_view(),
                                )
                            })
                        },
                        i,
                    )
                })
                .collect_view()
        },
    )
}

/// Poor man's lazy loading.
fn lazy<T>(view: impl FnOnce() -> T + Send + Sync + 'static, i: usize) -> impl IntoView
where
    T: IntoView + 'static,
{
    let trigger = ArcTrigger::new();
    let task = Task::new_local({
        let trigger = trigger.clone();
        async move {
            utile::time::sleep(Duration::from_millis((i as u64) * 10)).await;
            trigger.notify();
        }
    });
    let mut view = Some(view);
    let mut first = true;
    move || {
        trigger.track();
        let _task = &task;

        if first {
            first = false;
            None
        } else {
            Some(view.take().unwrap()())
        }
    }
}

fn render_edit_analysis_chart(
    edit_count: Signal<isize>,
    sample_count: usize,
    analysis: Signal<OutputGetEditAnalysis>,
) -> impl IntoView {
    let sample_count = sample_count as f64;
    let chart = Signal::derive_local(move || {
        // Untracked to avoid re-rendering before the data is ready. (Always trigger a change in the data.)
        let edit_count = edit_count.get_untracked();

        let counts: Vec<f64> = analysis
            .get()
            .edits
            .iter()
            .map(|variant_info| variant_info.count as f64)
            .map(|c| c / sample_count)
            .collect();

        let variant_count = analysis.get().edits.len();

        let mut plot = plotly::Plot::new();
        plot.add_trace(
            plotly::Histogram::new(counts.clone())
                .n_bins_x(100)
                .name("Edit Frequency Distribution"),
        );
        plot.set_layout(
            plotly::Layout::new()
                .template(PLOTLY_THEME.get())
                .x_axis(
                    plotly::layout::Axis::new()
                        .title(format!("Frequency in top {edit_count} edits"))
                        .tick_format(".0%"),
                )
                .y_axis(
                    plotly::layout::Axis::new()
                        .title(format!("Number of variants ({variant_count} total)")),
                )
                .margin(
                    plotly::layout::Margin::new()
                        .top(60)
                        .right(60)
                        .bottom(60)
                        .left(60),
                ),
        );
        plot
    });

    render_plotly(chart)
}

fn edit_table_line(
    variant_info: &VariantInfo,
    use_hq: bool,
    loading: ArcSignal<bool>,
) -> impl IntoView {
    let SummaryStats {
        chr,
        pos,
        ref_allele,
        alt,
        ..
    } = &*variant_info.stat;

    // Extract statistical information from the variant
    let VariantInfo {
        count: _,
        stat,
        ploidy_dosage,
    } = variant_info;

    // Format dosage/ploidy information
    let dosage_info = if ploidy_dosage.is_empty() {
        "No dosage data".to_owned()
    } else {
        ploidy_dosage
            .iter()
            .take(5) // Show top 5 most common dosage/ploidy combinations
            .map(|(VariantSampleInfo { dosage, ploidy }, count)| {
                let s = if *count == 1 { "" } else { "s" };
                format!("{dosage} of {ploidy} ({count} sample{s})")
            })
            .collect::<Vec<_>>()
            .join(", ")
    };

    let beta = if use_hq {
        let beta = stat.beta_meta_hq.unwrap();
        format!("β (HQ): {:.4}", beta)
    } else {
        let beta = stat.beta_meta.unwrap();
        format!("β: {:.4}", beta)
    };

    let (p_value_text, p_value) = if use_hq {
        let p_value = stat.neglog10_pval_meta_hq.unwrap();
        ("p-value (HQ):\u{00A0}10", -*p_value)
    } else {
        let p_value = stat.neglog10_pval_meta.unwrap();
        ("p-value:\u{00A0}10", -*p_value)
    };

    let gnomad_link = format!("https://gnomad.broadinstitute.org/variant/{chr}-{pos}-{ref_allele}-{alt}?dataset=gnomad_r4");

    html::div()
        .class("card edit-item")
        .class(("loading", loading))
        .child((
            html::div().class("edit-details").child((
                html::div().class("edit-primary").child((
                    html::span()
                        .class("edit-position")
                        .child(format!("{chr}:{pos}")),
                    html::span()
                        .class("edit-alleles")
                        .child(format!("{ref_allele} → {alt}")),
                    html::a()
                        .class("gnomad-link")
                        .attr("href", gnomad_link)
                        .attr("target", "_blank")
                        .attr("rel", "noopener noreferrer")
                        .child("gnomAD"),
                )),
                html::div().class("edit-stats").child(
                    html::span()
                        .class("stat-item")
                        .child(html::span().child((p_value_text, html::sup().child(p_value)))),
                ),
                html::div().class("edit-dosage").child((
                    html::span().class("dosage-label").child("Dosage×Ploidy: "),
                    html::span().class("dosage-info").child(dosage_info),
                )),
            )),
            html::div().child(html::span().class("count-badge").child(beta)),
        ))
}

fn render_plotly(
    plot: Signal<Plot, LocalStorage>,
) -> html::HtmlElement<html::Div, impl Attribute + Clone, impl RenderHtml + Clone> {
    let id = format!("plot_{}", rand::random::<u64>());

    Effect::new({
        let id = id.clone();
        move || {
            let mut plot = plot.get();
            plot.set_configuration(
                plot.configuration()
                    .clone()
                    .responsive(true)
                    .display_mode_bar(DisplayModeBar::False)
                    .display_logo(false)
                    .typeset_math(true),
            );

            let id = id.clone();
            let plot_update_task = async move {
                if document().get_element_by_id(&id).is_none() {
                    log::warn!("[App][Chart][{id}] Plot element not found");
                    return;
                }

                log::info!("[App][Chart][{id}] Updating plot");
                plotly::bindings::react(&id, &plot).await;
                log::info!("[App][Chart][{id}] Plot updated successfully");
            };
            // TODO: proper `defer` so it's cleaned up automatically.
            wasm_bindgen_futures::spawn_local(plot_update_task);
        }
    });

    html::div().id(id).style("width: 100%; height: 100%;")
}

fn toggle_button(
    signal: RwSignal<bool>,
    label: &'static str,
    tip: Option<&'static str>,
) -> impl IntoView {
    html::button()
        .class("toggle-btn")
        .class(("active", signal))
        .class(("tooltip", tip.is_some()))
        .attr("data-text", tip)
        .on(ev::click, move |_| signal.flip())
        .child(label)
}

impl PageState {
    fn new(file: String) -> Self {
        let initial_state = ArcRwSignal::new(None);

        let bridge = WorkerStruct::spawner()
            .callback({
                let initial_state = initial_state.clone();
                move |output| match output {
                    Output::Init(init) => initial_state.set(Some(init)),
                    Output::GetScores(_) | Output::GetEditAnalysis(_) => unreachable!(),
                }
            })
            .spawn("/worker.js");

        bridge.send(Input::Set {
            file: file.clone(),
            origin: leptos::prelude::window().location().origin().unwrap(),
        });

        let bridge = Arc::new(SendOption::new_local(Some(bridge)));

        Self {
            _file: file,

            bridge,

            init: initial_state,

            edit_count: RwSignal::new(0),
            use_hq: RwSignal::new(false),
            normalise: RwSignal::new(false),
            top_pvalues: RwSignal::new(100),

            show_stats: RwSignal::new(false),
            full_range: RwSignal::new(false),
        }
    }

    fn scores_signal(&self) -> Signal<Load<Option<ArcRwSignal<OutputGetScores>>>> {
        let bridge = self.bridge.clone();

        let edit_count = self.edit_count;
        let normalise = self.normalise;
        let use_hq_scores = self.use_hq;
        let top_pvalues = self.top_pvalues;

        let data = RwSignal::new(Load::Loading);

        Signal::derive(move || Input::GetScores {
            edit_count: edit_count.get(),
            normalise: normalise.get(),
            use_hq: use_hq_scores.get(),
            top_pvalues: top_pvalues.get(),
        })
        .rate_limit_leading(Duration::from_millis(10))
        .map_async(move |input| {
            let input = input.clone();
            let bridge = bridge.clone();
            async move {
                match query_bridge(&bridge, input).await {
                    Output::Init(_) | Output::GetEditAnalysis(_) => unreachable!(),
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

        data.into()
    }

    fn edits_signal(&self) -> Signal<Load<Option<ArcRwSignal<OutputGetEditAnalysis>>>> {
        let bridge = self.bridge.clone();

        let edit_count = self.edit_count;
        let normalise = self.normalise;
        let use_hq_scores = self.use_hq;
        let top_pvalues = self.top_pvalues;

        let data = RwSignal::new(Load::Loading);

        Signal::derive(move || Input::GetEditAnalysis {
            edit_count: edit_count.get(),
            normalise: normalise.get(),
            use_hq: use_hq_scores.get(),
            top_pvalues: top_pvalues.get(),
        })
        .rate_limit_leading(Duration::from_millis(10))
        .map_async(move |input| {
            let input = input.clone();
            let bridge = bridge.clone();
            async move {
                match query_bridge(&bridge, input).await {
                    Output::Init(_) | Output::GetScores(_) => unreachable!(),
                    Output::GetEditAnalysis(edits) => edits,
                }
            }
        })
        .for_each_immediate({
            let mut edits_signal: Option<ArcRwSignal<_>> = None;

            move |edits| match edits {
                Load::Loading => data.set(Load::Loading),
                Load::Ready(None) => data.set(Load::Ready(None)),
                Load::Ready(Some(edits)) => {
                    let edits = edits.clone();
                    let s = match edits_signal.clone() {
                        Some(s) => {
                            s.set(edits);
                            s
                        }
                        None => ArcRwSignal::new(edits),
                    };
                    edits_signal = Some(s.clone());
                    if !data.read_untracked().is_ready() {
                        data.set(Load::Ready(Some(s.clone())));
                    }
                }
            }
        });

        data.into()
    }
}

async fn query_bridge(bridge: &SendOption<WorkerBridge<WorkerStruct>>, input: Input) -> Output {
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
    let start = web_time::Instant::now();
    log::info!("[App] Requesting {input:?}");
    bridge.as_ref().unwrap().send(input.clone());
    let output = rx.await.unwrap();
    log::info!(
        "[App] Received output for {input:?} in {:?}ms",
        start.elapsed().as_millis()
    );
    output
}

fn document() -> web_sys::Document {
    web_sys::window().unwrap().document().unwrap()
}
