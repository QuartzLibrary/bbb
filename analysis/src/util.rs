use std::{cmp::Ordering, mem, ops::Deref, sync::Arc};

use biocore::{
    dna::{DnaBase, DnaSequence},
    location::{GenomePosition, SequenceOrientation},
    sequence::AsciiChar,
};
use hail::contig::GRCh38Contig;
use pan_ukbb::{PhenotypeManifestEntry, Population, SummaryStats};
use serde::{Deserialize, Serialize};

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

pub fn insert_ordered<T, F>(data: &mut Vec<T>, value: T, mut cmp: F)
where
    F: FnMut(&T, &T) -> Ordering,
{
    let i = data.partition_point(|v| cmp(v, &value).is_le());
    data.insert(i, value);
}

pub enum MArc<T> {
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
    pub fn new(t: T) -> Self {
        MArc::Inline(t)
    }
    pub fn arc(&mut self) -> Arc<T> {
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

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SummaryStatKey {
    pub chr: GRCh38Contig,
    pub pos: u64,
    pub ref_allele: DnaSequence,
    pub alt: DnaSequence,
}
impl std::fmt::Display for SummaryStatKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{} {}>{}",
            self.chr, self.pos, self.ref_allele, self.alt
        )
    }
}
impl std::str::FromStr for SummaryStatKey {
    type Err = std::io::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        fn invalid(e: impl Into<Box<dyn std::error::Error + Send + Sync>>) -> std::io::Error {
            std::io::Error::new(std::io::ErrorKind::InvalidInput, e)
        }
        let Some((chr_pos, ref_alt)) = s.split_once(' ') else {
            return Err(invalid("invalid summary stat key"));
        };
        let Some((chr, pos)) = chr_pos.split_once(':') else {
            return Err(invalid("invalid summary stat key"));
        };
        let Some((ref_allele, alt)) = ref_alt.split_once('>') else {
            return Err(invalid("invalid summary stat key"));
        };

        Ok(Self {
            chr: GRCh38Contig::from_str(chr).map_err(invalid)?,
            pos: pos.parse().map_err(invalid)?,
            ref_allele: DnaBase::decode(ref_allele.as_bytes().to_vec())?,
            alt: DnaBase::decode(alt.as_bytes().to_vec())?,
        })
    }
}
// We serialize as a string so that it can be used as a key in JSON.
impl Serialize for SummaryStatKey {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&self.to_string())
    }
}
impl<'de> Deserialize<'de> for SummaryStatKey {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        s.parse().map_err(serde::de::Error::custom)
    }
}
impl SummaryStatKey {
    pub fn new(summary_stat: &SummaryStats<GRCh38Contig>) -> Self {
        Self {
            chr: summary_stat.chr,
            pos: summary_stat.pos,
            ref_allele: summary_stat.ref_allele.clone(),
            alt: summary_stat.alt.clone(),
        }
    }

    pub fn at(&self) -> GenomePosition<GRCh38Contig> {
        GenomePosition {
            name: self.chr,
            orientation: SequenceOrientation::Forward,
            at: self.pos - 1,
        }
    }
}

pub trait InspectEvery: Iterator + Sized {
    fn inspect_every(
        self,
        n: usize,
        mut f: impl FnMut(usize, &<Self as Iterator>::Item),
    ) -> impl Iterator<Item = Self::Item> {
        let mut i = 0;
        self.inspect(move |item| {
            i += 1;
            if i % n == 0 {
                f(i, item);
            }
        })
    }
}
impl<I> InspectEvery for I where I: Iterator {}

pub struct FirstN<T, F> {
    n: usize,
    data: Vec<T>,
    cmp: F,
}
impl<T, F> Deref for FirstN<T, F> {
    type Target = [T];

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}
impl<T, F> FirstN<T, F>
where
    F: FnMut(&T, &T) -> Ordering,
{
    pub fn new(n: usize, cmp: F) -> Self {
        Self {
            n,
            data: Vec::new(),
            cmp,
        }
    }

    pub fn from_iter<I>(n: usize, cmp: F, iter: I) -> Self
    where
        I: IntoIterator<Item = T>,
    {
        let mut first_n = FirstN::new(n, cmp);
        for value in iter {
            first_n.push(value);
        }
        first_n
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn push(&mut self, value: T) {
        if self.n == 0 {
            return;
        }

        if self.data.len() < self.n {
            insert_ordered(&mut self.data, value, &mut self.cmp);
        } else {
            let last = self.data.last().unwrap();
            if (self.cmp)(&value, last).is_lt() {
                self.data.pop();
                insert_ordered(&mut self.data, value, &mut self.cmp);
            }
        }
    }

    pub fn push_ref(&mut self, value: &T)
    where
        T: Clone,
    {
        if self.n == 0 {
            return;
        }

        if self.data.len() < self.n {
            insert_ordered(&mut self.data, value.clone(), &mut self.cmp);
        } else {
            let last = self.data.last().unwrap();
            if (self.cmp)(value, last).is_lt() {
                self.data.pop();
                insert_ordered(&mut self.data, value.clone(), &mut self.cmp);
            }
        }
    }
}
impl<T, F> IntoIterator for FirstN<T, F> {
    type Item = T;
    type IntoIter = std::vec::IntoIter<T>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}
