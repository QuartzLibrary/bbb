use std::{
    cell::LazyCell,
    cmp::Ordering,
    collections::{BTreeSet, HashMap},
    ops::Range,
    sync::Arc,
};

use biocore::location::GenomePosition;
use genomes1000::{
    pedigree::{Pedigree, Sex},
    resource::Genomes1000Resource,
    simplified::SimplifiedRecord,
};
use hail::contig::GRCh38Contig;
use ordered_float::NotNan;
use pan_ukbb::{PhenotypeManifestEntry, SummaryStats};
use serde::{Deserialize, Serialize};
use std::hash::Hash;
use utile::resource::RawResourceExt;

use crate::util::{MArc, SummaryStatKey, mean, std_dev};

const EDIT_COUNT: usize = 100;

#[derive(Debug, Clone, PartialEq, Eq)]
#[derive(Serialize, Deserialize)]
pub struct Scores {
    pub phenotype: PhenotypeManifestEntry,
    pub scores: HashMap<String, Score>,
    /// Filled-in at the end before serialization.
    pub summary_stats: BTreeSet<SummaryStats<GRCh38Contig>>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[derive(Serialize, Deserialize)]
pub struct Score {
    pub score: NotNan<f64>,
    pub score_hq: NotNan<f64>,

    pub best: Vec<Variant>,
    pub worst: Vec<Variant>,
    pub best_hq: Vec<Variant>,
    pub worst_hq: Vec<Variant>,

    pub sex: Sex,
    pub population: String,
    pub superpopulation: String,
}
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[derive(Serialize, Deserialize)]
pub struct Variant {
    pub dosage: u8,
    pub ploidy: u8,

    pub max_edit: NotNan<f64>,
    pub min_edit: NotNan<f64>,

    pub max_edit_hq: Option<NotNan<f64>>,
    pub min_edit_hq: Option<NotNan<f64>>,

    pub key: Box<SummaryStatHandle>,
}

#[derive(Debug, Clone, Copy)]
#[derive(Serialize, Deserialize)]
pub struct Stats {
    pub mean: f64,
    pub std_dev: f64,
    pub min: f64,
    pub max: f64,
}
impl Stats {
    pub fn normalise_value(&self, v: f64) -> f64 {
        (v - self.mean) / self.std_dev
    }
    pub fn normalised_self(&self) -> Self {
        self.normalised(self.mean, self.std_dev)
    }
    pub fn normalised(&self, mean: f64, std_dev: f64) -> Self {
        let normalise = |v: f64| (v - mean) / std_dev;
        Self {
            mean: normalise(self.mean),
            std_dev: self.std_dev / std_dev,
            min: normalise(self.min),
            max: normalise(self.max),
        }
    }
}

impl Scores {
    pub async fn new(phenotype: PhenotypeManifestEntry) -> Self {
        let pedigrees = Genomes1000Resource::high_coverage_pedigree()
            .log_progress()
            .with_global_fs_cache()
            .ensure_cached_async()
            .await
            .unwrap();

        let mut scores: HashMap<String, Score> = HashMap::with_capacity(20_000);
        genomes1000::load_pedigree(pedigrees)
            .await
            .unwrap()
            .into_iter()
            .map(|p| (p.id.clone(), Score::new(p)))
            .collect_into(&mut scores);

        Self {
            phenotype,
            scores,
            summary_stats: BTreeSet::new(),
        }
    }

    pub fn push_missing(&mut self, summary_stat: SummaryStats<GRCh38Contig>) {
        let summary_stat = &mut MArc::new(summary_stat);
        for score in self.scores.values_mut() {
            score.push_missing(summary_stat);
        }
    }
    pub fn push_variant(
        &mut self,
        summary_stat: SummaryStats<GRCh38Contig>,
        sample_names: &[String],
        record: &SimplifiedRecord,
    ) {
        let summary_stat = &mut MArc::new(summary_stat);
        for (name, sample) in sample_names.iter().zip(&record.samples) {
            let score = self.scores.get_mut(name).unwrap();
            score.push_variant(
                summary_stat,
                sample.dosage(1),
                record.contig.ploidy(score.sex),
            );
        }
    }

    pub fn finalize(&mut self) {
        for score in self.scores.values_mut() {
            for variant in &mut score.best {
                self.summary_stats.insert(variant.summary_stat());
            }
            for variant in &mut score.worst {
                self.summary_stats.insert(variant.summary_stat());
            }
            for variant in &mut score.best_hq {
                self.summary_stats.insert(variant.summary_stat());
            }
            for variant in &mut score.worst_hq {
                self.summary_stats.insert(variant.summary_stat());
            }
        }
    }

    pub fn len(&self) -> usize {
        self.scores.len()
    }
    pub fn is_empty(&self) -> bool {
        self.scores.is_empty()
    }

    pub fn scores(&self) -> impl Iterator<Item = &Score> {
        self.scores.values()
    }

    pub fn full_score_range(&self) -> Option<Range<f64>> {
        let (min, max) = self.scores.values().map(|s| s.full_score_range()).fold(
            (None, None),
            |(min, max): (Option<f64>, Option<f64>), range| {
                let min = min.unwrap_or(range.start).min(range.start);
                let max = max.unwrap_or(range.end).max(range.end);
                (Some(min), Some(max))
            },
        );
        Some(min?..max?)
    }
    pub fn full_score_range_hq(&self) -> Option<Range<f64>> {
        let (min, max) = self.scores.values().map(|s| s.full_score_range_hq()).fold(
            (None, None),
            |(min, max): (Option<f64>, Option<f64>), range| {
                let min = min.unwrap_or(range.start).min(range.start);
                let max = max.unwrap_or(range.end).max(range.end);
                (Some(min), Some(max))
            },
        );
        Some(min?..max?)
    }

    pub fn edited_scores(&self, n: isize) -> impl Iterator<Item = Option<f64>> + use<'_> {
        self.scores.values().map(move |s| s.edited_score(n))
    }
    pub fn edited_scores_hq(&self, n: isize) -> impl Iterator<Item = Option<f64>> + use<'_> {
        self.scores.values().map(move |s| s.edited_score_hq(n))
    }

    pub fn max_edit_count(&self) -> isize {
        self.scores()
            .map(|s| s.best.len() as isize)
            .max()
            .unwrap_or(0)
    }
    pub fn min_edit_count(&self) -> isize {
        self.scores()
            .map(|s| s.worst.len() as isize)
            .min()
            .unwrap_or(0)
    }

    pub fn max_edit_count_hq(&self) -> isize {
        self.scores()
            .map(|s| s.best_hq.len() as isize)
            .max()
            .unwrap_or(0)
    }
    pub fn min_edit_count_hq(&self) -> isize {
        self.scores()
            .map(|s| s.worst_hq.len() as isize)
            .min()
            .unwrap_or(0)
    }

    pub fn stats(&self) -> Option<Stats> {
        Some(Stats {
            mean: *self.mean()?,
            std_dev: *self.std_dev().unwrap(),
            min: *self.min().unwrap(),
            max: *self.max().unwrap(),
        })
    }
    pub fn mean(&self) -> Option<NotNan<f64>> {
        Some(NotNan::new(mean(self.scores().map(|s| *s.score))?).unwrap())
    }
    pub fn std_dev(&self) -> Option<NotNan<f64>> {
        Some(NotNan::new(std_dev(*self.mean()?, self.scores().map(|s| *s.score))?).unwrap())
    }
    pub fn min(&self) -> Option<NotNan<f64>> {
        self.scores().map(|s| s.score).min()
    }
    pub fn max(&self) -> Option<NotNan<f64>> {
        self.scores().map(|s| s.score).max()
    }

    pub fn stats_hq(&self) -> Option<Stats> {
        Some(Stats {
            mean: *self.mean_hq()?,
            std_dev: *self.std_dev_hq().unwrap(),
            min: *self.min_hq().unwrap(),
            max: *self.max_hq().unwrap(),
        })
    }
    pub fn mean_hq(&self) -> Option<NotNan<f64>> {
        Some(NotNan::new(mean(self.scores().map(|s| *s.score_hq))?).unwrap())
    }
    pub fn std_dev_hq(&self) -> Option<NotNan<f64>> {
        Some(
            NotNan::new(std_dev(
                *self.mean_hq()?,
                self.scores().map(|s| *s.score_hq),
            )?)
            .unwrap(),
        )
    }
    pub fn min_hq(&self) -> Option<NotNan<f64>> {
        self.scores().map(|s| s.score_hq).min()
    }
    pub fn max_hq(&self) -> Option<NotNan<f64>> {
        self.scores().map(|s| s.score_hq).max()
    }

    pub fn normalised(&self) -> Self {
        if self.is_empty() {
            return self.clone();
        }

        let mean = self.mean().unwrap();
        let std_dev = self.std_dev().unwrap();

        let mean_hq = self.mean_hq().unwrap();
        let std_dev_hq = self.std_dev_hq().unwrap();

        Self {
            phenotype: self.phenotype.clone(),
            scores: self
                .scores
                .clone()
                .into_iter()
                .map(|(k, v)| (k, v.normalised(mean, std_dev, mean_hq, std_dev_hq)))
                .collect(),
            summary_stats: self
                .summary_stats
                .clone()
                .into_iter()
                .map(|s| s.normalised(std_dev, std_dev_hq))
                .collect(),
        }
    }
}

impl Score {
    fn new(pedigree: Pedigree) -> Self {
        Self {
            score: NotNan::new(0.0).unwrap(),
            score_hq: NotNan::new(0.0).unwrap(),

            best: Vec::with_capacity(EDIT_COUNT + 1),
            worst: Vec::with_capacity(EDIT_COUNT + 1),

            best_hq: Vec::with_capacity(EDIT_COUNT + 1),
            worst_hq: Vec::with_capacity(EDIT_COUNT + 1),

            sex: pedigree.sex,
            population: pedigree.population,
            superpopulation: pedigree.superpopulation,
        }
    }
    fn push_missing(&mut self, s: &mut MArc<SummaryStats<GRCh38Contig>>) {
        let ploidy = match (self.sex, s.chr) {
            (Sex::Male, GRCh38Contig::X) => 1,
            (Sex::Female, GRCh38Contig::X) => 2,
            (Sex::Male, GRCh38Contig::Y) => 1,
            (Sex::Female, GRCh38Contig::Y) => 0,
            (_, GRCh38Contig::MT) => 1,
            (_, _) => 2,
        };
        self.push_variant(s, 0, ploidy);
    }
    fn push_variant(&mut self, s: &mut MArc<SummaryStats<GRCh38Contig>>, dosage: u8, ploidy: u8) {
        {
            let s = &**s;
            assert!(
                dosage <= ploidy,
                "dosage: {dosage}, ploidy: {ploidy}, s: {s:?}"
            );
            assert!(ploidy <= 2, "dosage: {dosage}, ploidy: {ploidy}, s: {s:?}");
        }
        self.score += s.beta_meta.unwrap() * NotNan::from(dosage);
        if let Some(beta_meta_hq) = s.beta_meta_hq {
            self.score_hq += beta_meta_hq * NotNan::from(dosage);
        }

        let max = s.max_edit(dosage, ploidy).unwrap();
        let min = s.min_edit(dosage, ploidy).unwrap();
        let max_hq = s.max_edit_hq(dosage, ploidy);
        let min_hq = s.min_edit_hq(dosage, ploidy);

        let v = LazyCell::new(move || Variant::new(s.arc(), dosage, ploidy));

        if self.best_bound() < max {
            let i = self.best.partition_point(|v| v.max_edit >= max);
            if i < EDIT_COUNT {
                self.best.insert(i, v.clone());
                self.best.truncate(EDIT_COUNT);
            }
        }
        if min < self.worst_bound() {
            let i = self.worst.partition_point(|v| v.min_edit <= min);
            if i < EDIT_COUNT {
                self.worst.insert(i, v.clone());
                self.worst.truncate(EDIT_COUNT);
            }
        }
        if let Some(max_hq) = max_hq {
            if self.best_bound_hq() < max_hq {
                let i = self
                    .best_hq
                    .partition_point(|v| v.max_edit_hq.unwrap() >= max_hq);
                if i < EDIT_COUNT {
                    self.best_hq.insert(i, v.clone());
                    self.best_hq.truncate(EDIT_COUNT);
                }
            }
        }
        if let Some(min_hq) = min_hq {
            if min_hq < self.worst_bound_hq() {
                let i = self
                    .worst_hq
                    .partition_point(|v| v.min_edit_hq.unwrap() <= min_hq);
                if i < EDIT_COUNT {
                    self.worst_hq.insert(i, v.clone());
                    self.worst_hq.truncate(EDIT_COUNT);
                }
            }
        }
    }
    fn best_bound(&self) -> NotNan<f64> {
        let zero = NotNan::new(0.).unwrap();
        if self.best.len() < EDIT_COUNT {
            zero
        } else {
            self.best.last().map(|v| v.max_edit).unwrap_or(zero)
        }
    }
    fn worst_bound(&self) -> NotNan<f64> {
        let zero = NotNan::new(0.).unwrap();
        if self.worst.len() < EDIT_COUNT {
            zero
        } else {
            self.worst.last().map(|v| v.min_edit).unwrap_or(zero)
        }
    }

    fn best_bound_hq(&self) -> NotNan<f64> {
        let zero = NotNan::new(0.).unwrap();
        if self.best_hq.len() < EDIT_COUNT {
            zero
        } else {
            self.best_hq
                .last()
                .map(|v| v.max_edit_hq.unwrap_or(zero))
                .unwrap_or(zero)
        }
    }
    fn worst_bound_hq(&self) -> NotNan<f64> {
        let zero = NotNan::new(0.).unwrap();
        if self.worst_hq.len() < EDIT_COUNT {
            zero
        } else {
            self.worst_hq
                .last()
                .map(|v| v.min_edit_hq.unwrap_or(zero))
                .unwrap_or(zero)
        }
    }
}
impl Score {
    pub fn full_score_range(&self) -> Range<f64> {
        let min = self.edited_score(-(self.worst.len() as isize)).unwrap();
        let max = self.edited_score(self.best.len() as isize).unwrap();
        min..max
    }
    pub fn full_score_range_hq(&self) -> Range<f64> {
        let min = self
            .edited_score_hq(-(self.worst_hq.len() as isize))
            .unwrap();
        let max = self.edited_score_hq(self.best_hq.len() as isize).unwrap();
        min..max
    }

    pub fn edited_score(&self, n: isize) -> Option<f64> {
        Some(
            self.score
                + match Ord::cmp(&n, &0) {
                    Ordering::Greater => self
                        .best
                        .get(..n as usize)?
                        .iter()
                        .map(|v| *v.max_edit)
                        .sum::<f64>(),
                    Ordering::Less => self
                        .worst
                        .get(..-n as usize)?
                        .iter()
                        .map(|v| *v.min_edit)
                        .sum::<f64>(),
                    Ordering::Equal => 0.,
                },
        )
    }
    pub fn edited_score_hq(&self, n: isize) -> Option<f64> {
        Some(
            self.score_hq
                + match Ord::cmp(&n, &0) {
                    Ordering::Greater => self
                        .best_hq
                        .get(..n as usize)?
                        .iter()
                        .map(|v| *v.max_edit_hq.unwrap())
                        .sum::<f64>(),
                    Ordering::Less => self
                        .worst_hq
                        .get(..-n as usize)?
                        .iter()
                        .map(|v| *v.min_edit_hq.unwrap())
                        .sum::<f64>(),
                    Ordering::Equal => 0.,
                },
        )
    }

    pub fn normalised(
        &self,
        mean: NotNan<f64>,
        std_dev: NotNan<f64>,
        mean_hq: NotNan<f64>,
        std_dev_hq: NotNan<f64>,
    ) -> Self {
        Self {
            score: (self.score - mean) / std_dev,
            score_hq: (self.score_hq - mean_hq) / std_dev_hq,

            best: self
                .best
                .iter()
                .map(|v| v.normalised(std_dev, std_dev_hq))
                .collect(),
            worst: self
                .worst
                .iter()
                .map(|v| v.normalised(std_dev, std_dev_hq))
                .collect(),

            best_hq: self
                .best_hq
                .iter()
                .map(|v| v.normalised(std_dev, std_dev_hq))
                .collect(),
            worst_hq: self
                .worst_hq
                .iter()
                .map(|v| v.normalised(std_dev, std_dev_hq))
                .collect(),

            ..self.clone()
        }
    }
}
impl Variant {
    pub fn new(s: Arc<SummaryStats<GRCh38Contig>>, dosage: u8, ploidy: u8) -> Self {
        Self {
            dosage,
            ploidy,

            max_edit: s.max_edit(dosage, ploidy).unwrap(),
            min_edit: s.min_edit(dosage, ploidy).unwrap(),

            max_edit_hq: s.max_edit_hq(dosage, ploidy),
            min_edit_hq: s.min_edit_hq(dosage, ploidy),

            key: Box::new(SummaryStatHandle::new(s)),
        }
    }
    pub fn summary_stat(&self) -> SummaryStats<GRCh38Contig> {
        let s: &SummaryStats<GRCh38Contig> = self.key._handle.as_ref().unwrap();
        s.clone()
    }

    pub fn normalised(&self, std_dev: NotNan<f64>, std_dev_hq: NotNan<f64>) -> Self {
        Self {
            max_edit: self.max_edit / std_dev,
            min_edit: self.min_edit / std_dev,

            max_edit_hq: self.max_edit_hq.map(|v| v / std_dev_hq),
            min_edit_hq: self.min_edit_hq.map(|v| v / std_dev_hq),

            ..self.clone()
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[derive(Deserialize, Serialize)]
pub struct SummaryStatHandle {
    #[serde(flatten, with = "legacy")]
    pub key: SummaryStatKey,
    #[serde(skip)]
    pub _handle: Option<Arc<SummaryStats<GRCh38Contig>>>,
}
impl SummaryStatHandle {
    pub fn new(summary_stat: Arc<SummaryStats<GRCh38Contig>>) -> Self {
        Self {
            key: SummaryStatKey::new(&summary_stat),
            _handle: Some(summary_stat),
        }
    }
    pub fn new_handleless(summary_stat: &SummaryStats<GRCh38Contig>) -> Self {
        Self {
            key: SummaryStatKey::new(summary_stat),
            _handle: None,
        }
    }

    pub fn at(&self) -> GenomePosition<GRCh38Contig> {
        self.key.at()
    }
}

mod legacy {
    use biocore::dna::DnaSequence;
    use hail::contig::GRCh38Contig;
    use serde::{Deserialize, Serialize};

    use crate::util::SummaryStatKey;

    #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
    #[derive(Deserialize)]
    pub struct SummaryStatKeyCopy {
        pub chr: GRCh38Contig,
        pub pos: u64,
        pub ref_allele: DnaSequence,
        pub alt: DnaSequence,
    }

    #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
    #[derive(Serialize)]
    pub struct SummaryStatKeyRef<'a> {
        pub chr: GRCh38Contig,
        pub pos: u64,
        pub ref_allele: &'a DnaSequence,
        pub alt: &'a DnaSequence,
    }

    impl<'a> From<&'a SummaryStatKey> for SummaryStatKeyRef<'a> {
        fn from(key: &'a SummaryStatKey) -> Self {
            Self {
                chr: key.chr,
                pos: key.pos,
                ref_allele: &key.ref_allele,
                alt: &key.alt,
            }
        }
    }
    impl From<SummaryStatKeyCopy> for SummaryStatKey {
        fn from(key: SummaryStatKeyCopy) -> Self {
            Self {
                chr: key.chr,
                pos: key.pos,
                ref_allele: key.ref_allele,
                alt: key.alt,
            }
        }
    }

    pub fn serialize<S>(v: &SummaryStatKey, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let key = SummaryStatKeyRef::from(v);
        key.serialize(serializer)
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<SummaryStatKey, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let key: SummaryStatKeyCopy = Deserialize::deserialize(deserializer)?;
        Ok(key.into())
    }
}
