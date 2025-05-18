use pan_ukbb::{PhenotypeManifestEntry, Population};

pub fn mean(value: impl Iterator<Item = f64>) -> Option<f64> {
    let mut sum = 0.0;
    let mut count = 0.0;
    for v in value {
        sum += v;
        count += 1.0;
    }
    if count > 0. { Some(sum / count) } else { None }
}

pub fn std_dev(mean: f64, value: impl Iterator<Item = f64>) -> Option<f64> {
    let mut sum_of_squares = 0.0;
    let mut count = 0.0;
    for v in value {
        sum_of_squares += (v - mean).powi(2);
        count += 1.0;
    }
    if count > 0. {
        Some((sum_of_squares / count).sqrt())
    } else {
        None
    }
}

pub fn passes_qc(phenotype: &PhenotypeManifestEntry) -> bool {
    // Kind of arbitrary, just to shave the number.
    phenotype.pops_pass_qc.contains(&Population::Eur)
        && phenotype.n_cases_full_cohort_both_sexes > 50_000
}
