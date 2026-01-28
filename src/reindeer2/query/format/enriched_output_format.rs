use super::super::load_kmer_counts_vector;
use crate::reindeer2::{BreakpointsNormalize, MatrixFormat, OutputFormat};

#[derive(Clone, Debug, PartialEq)]
pub struct KmerCountsAndNormalizeValue {
    pub kmer_counts: Vec<usize>,
    pub normalize_value: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BreakpointsXorEnrichedNormalize {
    Breakpoints(f64),
    Normalize(KmerCountsAndNormalizeValue),
}

impl BreakpointsNormalize {
    fn enriched(&self, bf_dir: &str) -> BreakpointsXorEnrichedNormalize {
        match self {
            BreakpointsNormalize::Breakpoints(penalty) => {
                BreakpointsXorEnrichedNormalize::Breakpoints(*penalty)
            }
            BreakpointsNormalize::Normalize(normalized) => {
                BreakpointsXorEnrichedNormalize::Normalize(load_kmer_count(*normalized, bf_dir))
            }
        }
    }
}

// TODO rename
#[derive(Clone, Debug, PartialEq)]
pub enum EnrichedMatrixFormat {
    Raw(
        /// None if not computing breakpoints nor normalize
        Option<BreakpointsXorEnrichedNormalize>,
    ),
    Average {
        normalized: Option<KmerCountsAndNormalizeValue>,
    },
    Median {
        normalized: Option<KmerCountsAndNormalizeValue>,
    },
}

impl EnrichedMatrixFormat {
    fn add_kmer_count(format_infos: MatrixFormat, bf_dir: &str) -> Self {
        match format_infos {
            MatrixFormat::Raw(raw_format_infos) => {
                Self::Raw(raw_format_infos.map(|x| x.enriched(bf_dir)))
            }
            MatrixFormat::Average { normalized } => Self::Average {
                normalized: load_kmer_count_from_option(normalized, bf_dir),
            },
            MatrixFormat::Median { normalized } => Self::Median {
                normalized: load_kmer_count_from_option(normalized, bf_dir),
            },
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum EnrichedOutputFormat {
    Colored {
        normalized: Option<KmerCountsAndNormalizeValue>,
    },
    Median {
        normalized: Option<KmerCountsAndNormalizeValue>,
    },
    AbundanceMatrix {
        format: EnrichedMatrixFormat,
    },
}

// TODO rename ?
fn load_kmer_count(normalize_value: u64, bf_dir: &str) -> KmerCountsAndNormalizeValue {
    KmerCountsAndNormalizeValue {
        kmer_counts: load_kmer_counts_vector(bf_dir)
            .expect("Failed to load from disk the kmer counts vector"),
        normalize_value,
    }
}

fn load_kmer_count_from_option(
    normalized: Option<u64>,
    bf_dir: &str,
) -> Option<KmerCountsAndNormalizeValue> {
    // we need the count of kmers if we want to normalize them
    normalized.map(|normalize_value| load_kmer_count(normalize_value, bf_dir))
}

impl EnrichedOutputFormat {
    pub fn from_pub_output_format(pub_output_format: OutputFormat, bf_dir: &str) -> Self {
        match pub_output_format {
            OutputFormat::AbundanceMatrix { format } => Self::AbundanceMatrix {
                format: EnrichedMatrixFormat::add_kmer_count(format, bf_dir),
            },
            OutputFormat::Colored { normalized } => Self::Colored {
                normalized: load_kmer_count_from_option(normalized, bf_dir),
            },
            OutputFormat::Median { normalized } => Self::Median {
                normalized: load_kmer_count_from_option(normalized, bf_dir),
            },
        }
    }
}
