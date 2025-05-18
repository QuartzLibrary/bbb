use std::collections::BTreeMap;

use genomes1000::{
    pedigree::{Pedigree, Sex},
    resource::Genomes1000Resource,
    simplified::SimplifiedRecord,
};
use hail::contig::GRCh38Contig;
use ordered_float::NotNan;
use pan_ukbb::{PhenotypeManifestEntry, SummaryStats};
use serde::{Deserialize, Serialize};
use utile::resource::RawResourceExt;

use crate::util::SummaryStatKey;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[derive(Serialize, Deserialize)]
pub struct PhenotypeTopPValues {
    pub phenotype: PhenotypeManifestEntry,
    pub top_variants: BTreeMap<SummaryStatKey, TopPValueVariant>,
    pub samples: BTreeMap<String, Pedigree>,
}
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[derive(Serialize, Deserialize)]
pub struct TopPValueVariant {
    pub stat: SummaryStats<GRCh38Contig>,

    pub genotypes: BTreeMap<String, SampleGenotype>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[derive(Serialize, Deserialize)]
pub struct SampleGenotype {
    pub dosage: u8,
    pub ploidy: u8,
}

impl PhenotypeTopPValues {
    pub async fn new(phenotype: PhenotypeManifestEntry) -> Self {
        let pedigrees = Genomes1000Resource::high_coverage_pedigree()
            .log_progress()
            .with_global_fs_cache()
            .ensure_cached_async()
            .await
            .unwrap();

        let samples: BTreeMap<String, Pedigree> = genomes1000::load_pedigree(pedigrees)
            .await
            .unwrap()
            .into_iter()
            .map(|p| (p.id.clone(), p))
            .collect();

        Self {
            phenotype,
            top_variants: BTreeMap::new(),
            samples,
        }
    }

    pub fn push_missing(&mut self, summary_stat: SummaryStats<GRCh38Contig>) {
        self.top_variants.insert(
            SummaryStatKey::new(&summary_stat),
            TopPValueVariant::new_missing(
                summary_stat,
                self.samples.iter().map(|(name, p)| (name.clone(), p.sex)),
            ),
        );
    }
    pub fn push_variant(
        &mut self,
        summary_stat: SummaryStats<GRCh38Contig>,
        sample_names: &[String],
        record: &SimplifiedRecord,
    ) {
        self.top_variants.insert(
            SummaryStatKey::new(&summary_stat),
            TopPValueVariant::new(
                summary_stat,
                sample_names,
                |name| self.samples[name].sex,
                record,
            ),
        );
    }

    pub fn len(&self) -> usize {
        self.top_variants.len()
    }
    pub fn is_empty(&self) -> bool {
        self.top_variants.is_empty()
    }

    pub fn normalised(&self, std_dev: NotNan<f64>, std_dev_hq: NotNan<f64>) -> Self {
        if self.is_empty() {
            return self.clone();
        }

        Self {
            phenotype: self.phenotype.clone(),
            top_variants: self
                .top_variants
                .clone()
                .into_iter()
                .map(|(k, v)| (k, v.normalised(std_dev, std_dev_hq)))
                .collect(),
            samples: self.samples.clone(),
        }
    }
}

impl TopPValueVariant {
    fn new_missing(
        summary_stat: SummaryStats<GRCh38Contig>,
        samples: impl Iterator<Item = (String, Sex)>,
    ) -> Self {
        let genotypes = samples
            .map(|(name, sex)| {
                (
                    name,
                    SampleGenotype {
                        dosage: 0,
                        ploidy: summary_stat.chr.ploidy(sex == Sex::Male),
                    },
                )
            })
            .collect();
        Self {
            stat: summary_stat,
            genotypes,
        }
    }
    fn new(
        summary_stat: SummaryStats<GRCh38Contig>,
        sample_names: &[String],
        get_sex: impl Fn(&str) -> Sex,
        record: &SimplifiedRecord,
    ) -> Self {
        let genotypes = sample_names
            .iter()
            .cloned()
            .zip(&record.samples)
            .map(|(name, sample)| {
                let sex = get_sex(&name);
                (
                    name,
                    SampleGenotype {
                        dosage: sample.dosage(1),
                        ploidy: summary_stat.chr.ploidy(sex == Sex::Male),
                    },
                )
            })
            .collect();
        Self {
            stat: summary_stat,
            genotypes,
        }
    }
}
impl TopPValueVariant {
    pub fn normalised(&self, std_dev: NotNan<f64>, std_dev_hq: NotNan<f64>) -> Self {
        Self {
            stat: self.stat.clone().normalised(std_dev, std_dev_hq),
            genotypes: self.genotypes.clone(),
        }
    }
}
