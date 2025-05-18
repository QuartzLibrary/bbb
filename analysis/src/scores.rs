use std::{
    cell::LazyCell,
    cmp::Ordering,
    collections::{BTreeMap, BTreeSet},
    mem,
    ops::{Deref, Range},
    sync::Arc,
};

use biocore::dna::DnaSequence;
use genomes1000::{
    pedigree::{Pedigree, Sex},
    resource::Genomes1000Resource,
    simplified::SimplifiedRecord,
};
use ordered_float::NotNan;
use pan_ukbb::{PhenotypeManifestEntry, SummaryStats};
use serde::{Deserialize, Serialize};
use std::hash::Hash;
use utile::resource::RawResourceExt;

use crate::util::{mean, std_dev};

const EDIT_COUNT: usize = 500;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[derive(Serialize, Deserialize)]
pub struct Scores {
    pub phenotype: PhenotypeManifestEntry,
    pub scores: BTreeMap<String, Score>,
    /// Filled-in at the end before serialization.
    pub summary_stats: BTreeSet<SummaryStats>,
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

    pub key: SummaryStatKey,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[derive(Serialize, Deserialize)]
pub struct SummaryStatKey {
    pub chr: String,
    pub pos: u64,
    pub ref_allele: DnaSequence,
    pub alt: DnaSequence,
    #[serde(skip)]
    pub _handle: Option<Arc<SummaryStats>>,
}

#[derive(Debug, Clone, Copy)]
#[derive(Serialize, Deserialize)]
pub struct Stats {
    pub mean: f64,
    pub std_dev: f64,
    pub min: f64,
    pub max: f64,
}

impl Scores {
    pub async fn new(phenotype: PhenotypeManifestEntry) -> Self {
        let pedigrees = Genomes1000Resource::high_coverage_pedigree()
            .log_progress()
            .with_global_fs_cache()
            .ensure_cached_async()
            .await
            .unwrap();

        let scores: BTreeMap<String, Score> = genomes1000::load_pedigree(pedigrees)
            .await
            .unwrap()
            .into_iter()
            .map(|p| (p.id.clone(), Score::new(p)))
            .collect();

        Self {
            phenotype,
            scores,
            summary_stats: BTreeSet::new(),
        }
    }

    pub fn push_missing(&mut self, summary_stat: SummaryStats) {
        let summary_stat = &mut MArc::new(summary_stat);
        for score in self.scores.values_mut() {
            score.push_missing(summary_stat);
        }
    }
    pub fn push_variant(
        &mut self,
        summary_stat: SummaryStats,
        sample_names: &[String],
        record: &SimplifiedRecord,
    ) {
        let summary_stat = &mut MArc::new(summary_stat);
        for (name, sample) in sample_names.iter().zip(&record.samples) {
            let score = self.scores.get_mut(name).unwrap();
            score.push_variant(
                summary_stat,
                sample.dosage(),
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

    pub fn full_score_range(&self) -> Range<f64> {
        let (min, max) = self.scores.values().map(|s| s.full_score_range()).fold(
            (None, None),
            |(min, max): (Option<f64>, Option<f64>), range| {
                let min = min.unwrap_or(range.start).min(range.start);
                let max = max.unwrap_or(range.end).max(range.end);
                (Some(min), Some(max))
            },
        );
        min.unwrap_or(0.)..max.unwrap_or(0.)
    }
    pub fn full_score_range_hq(&self) -> Range<f64> {
        let (min, max) = self.scores.values().map(|s| s.full_score_range_hq()).fold(
            (None, None),
            |(min, max): (Option<f64>, Option<f64>), range| {
                let min = min.unwrap_or(range.start).min(range.start);
                let max = max.unwrap_or(range.end).max(range.end);
                (Some(min), Some(max))
            },
        );
        min.unwrap_or(0.)..max.unwrap_or(0.)
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
                .map(|s| s.normalised(mean, std_dev, mean_hq, std_dev_hq))
                .collect(),
        }
    }
}

impl Score {
    fn new(pedigree: Pedigree) -> Self {
        Self {
            score: NotNan::new(0.0).unwrap(),
            score_hq: NotNan::new(0.0).unwrap(),

            best: Vec::new(),
            worst: Vec::new(),

            best_hq: Vec::new(),
            worst_hq: Vec::new(),

            sex: pedigree.sex,
            population: pedigree.population,
            superpopulation: pedigree.superpopulation,
        }
    }
    fn push_missing(&mut self, s: &mut MArc<SummaryStats>) {
        let ploidy = match (self.sex, &*s.chr) {
            (Sex::Male, "chrX") => 1,
            (Sex::Female, "chrX") => 2,
            (Sex::Male, "chrY") => 1,
            (Sex::Female, "chrY") => 0,
            (_, "chrM") => 1,
            (_, _) => 2,
        };
        self.push_variant(s, 0, ploidy);
    }
    fn push_variant(&mut self, s: &mut MArc<SummaryStats>, dosage: u8, ploidy: u8) {
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
        let max = NotNan::new(f64::MAX).unwrap();
        if self.best_hq.len() < EDIT_COUNT {
            max
        } else {
            self.best_hq
                .last()
                .map(|v| v.max_edit_hq.unwrap_or(max))
                .unwrap_or(max)
        }
    }
    fn worst_bound_hq(&self) -> NotNan<f64> {
        let min = NotNan::new(f64::MIN).unwrap();
        if self.worst_hq.len() < EDIT_COUNT {
            min
        } else {
            self.worst_hq
                .last()
                .map(|v| v.min_edit_hq.unwrap_or(min))
                .unwrap_or(min)
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
                .map(|v| v.normalised(std_dev_hq, std_dev_hq))
                .collect(),
            worst_hq: self
                .worst_hq
                .iter()
                .map(|v| v.normalised(std_dev_hq, std_dev_hq))
                .collect(),

            ..self.clone()
        }
    }
}
impl Variant {
    fn new(s: Arc<SummaryStats>, dosage: u8, ploidy: u8) -> Self {
        let beta_meta = s.beta_meta.unwrap();

        let edit_force_alt = NotNan::from(ploidy - dosage) * beta_meta;
        let edit_force_ref = -NotNan::from(dosage) * beta_meta;

        let edit_force_alt_hq = s
            .beta_meta_hq
            .map(|beta_meta_hq| NotNan::from(ploidy - dosage) * beta_meta_hq);
        let edit_force_ref_hq = s
            .beta_meta_hq
            .map(|beta_meta_hq| -NotNan::from(dosage) * beta_meta_hq);

        Self {
            dosage,
            ploidy,

            max_edit: Ord::max(edit_force_alt, edit_force_ref),
            min_edit: Ord::min(edit_force_alt, edit_force_ref),

            max_edit_hq: Ord::max(edit_force_alt_hq, edit_force_ref_hq),
            min_edit_hq: Ord::min(edit_force_alt_hq, edit_force_ref_hq),

            key: SummaryStatKey::new(s),
        }
    }
    pub fn summary_stat(&self) -> SummaryStats {
        let s: &SummaryStats = self.key._handle.as_ref().unwrap();
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

impl SummaryStatKey {
    pub fn new(summary_stat: Arc<SummaryStats>) -> Self {
        Self {
            chr: summary_stat.chr.clone(),
            pos: summary_stat.pos,
            ref_allele: summary_stat.ref_allele.clone(),
            alt: summary_stat.alt.clone(),
            _handle: Some(summary_stat),
        }
    }
}

enum MArc<T> {
    Invalid,
    Inline(T),
    Arc(Arc<T>),
}
impl<T> Deref for MArc<T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        match self {
            MArc::Invalid => unreachable!("Invalid"),
            MArc::Inline(t) => t,
            MArc::Arc(t) => t,
        }
    }
}
impl<T> MArc<T> {
    fn new(t: T) -> Self {
        MArc::Inline(t)
    }
    fn arc(&mut self) -> Arc<T> {
        match mem::replace(self, MArc::Invalid) {
            MArc::Invalid => unreachable!("Invalid"),
            MArc::Inline(t) => {
                let arc = Arc::new(t);
                *self = MArc::Arc(arc.clone());
                arc
            }
            MArc::Arc(t) => {
                *self = MArc::Arc(t.clone());
                t
            }
        }
    }
}
