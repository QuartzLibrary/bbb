use leptos::{
    html::{self},
    prelude::*,
    IntoView,
};

use pan_ukbb::PhenotypeManifestEntry;

pub fn home() -> impl IntoView {
    html::div()
        .class("manifest")
        .child(html::h1().child("Manifest"))
        .child(|| {
            let manifest: ArcRwSignal<_> = crate::MANIFEST.with(|m| (**m).clone());
            manifest.with(|manifest| match manifest {
                Some(Ok(manifest)) => html::div()
                    .child(
                        manifest
                            .iter()
                            .filter(|p| analysis::util::passes_qc(p))
                            .map(manifest_entry)
                            .collect_view(),
                    )
                    .into_any(),
                Some(Err(e)) => {
                    log::error!("[App] Error loading manifest: {}", e);
                    html::div()
                        .child("Error loading the manifest...")
                        .into_any()
                }
                None => html::div().child("Loading...").into_any(),
            })
        })
}
fn manifest_entry(entry: &PhenotypeManifestEntry) -> impl IntoView {
    let PhenotypeManifestEntry {
        trait_type,
        phenocode,
        pheno_sex,
        description,
        description_more,
        category,
        n_cases_full_cohort_both_sexes,
        n_cases_hq_cohort_both_sexes,
        filename,
        ..
    } = entry;

    let identifier = filename.strip_suffix(".tsv.bgz").unwrap_or(filename);
    let url = format!("/pan_ukbb/{identifier}");

    html::div()
        .class("manifest-entry")
        .child(
            html::div().class("entry-header").child(
                html::a()
                    .class("phenotype-title")
                    .href(url)
                    .child(description.clone()),
            ),
        )
        .child(
            html::div().class("entry-details").child((
                html::div().class("detail-row").child((
                    html::a()
                        .href(format!(
                            "https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id={phenocode}"
                        ))
                        .child(format!("[{phenocode}]")),
                    format!(" ({trait_type}) ({pheno_sex})"),
                )),
                html::div().class("detail-row").child((
                    html::span().class("detail-label").child("Category: "),
                    html::span()
                        .class("detail-value")
                        .child(category.clone().unwrap_or_default()),
                )),
                html::div().class("detail-row").child((
                    html::span().class("detail-label").child("Cases: "),
                    html::span().class("detail-value").child(format!(
                        "{n_cases_full_cohort_both_sexes} ({})",
                        n_cases_hq_cohort_both_sexes.unwrap_or(0)
                    )),
                )),
                description_more.clone().map(|d| {
                    html::div().class("detail-row").child((
                        html::span().class("detail-label").child("Description: "),
                        html::span().class("detail-value").child(d),
                    ))
                }),
            )),
        )
}
