use super::super::load_kmer_counts_vector;
use crate::reindeer2::OutputFormat;

#[derive(Clone, Debug, PartialEq)]
pub struct KmerCountsAndNormalizeValue {
    pub kmer_counts: Vec<usize>,
    pub normalize_value: u64,
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
        normalized: Option<KmerCountsAndNormalizeValue>,

        breakpoints: Option<f64>,
    },
}

impl EnrichedOutputFormat {
    pub fn add_kmer_count(
        normalized: Option<u64>,
        bf_dir: &str,
    ) -> std::option::Option<KmerCountsAndNormalizeValue> {
        // we need the count of kmers if we want to normalize them
        normalized.map(|normalize_value| KmerCountsAndNormalizeValue {
            kmer_counts: load_kmer_counts_vector(bf_dir)
                .expect("Failed to load from disk the kmer counts vector"),
            normalize_value,
        })
    }

    pub fn from_pub_output_format(pub_output_format: OutputFormat, bf_dir: &str) -> Self {
        match pub_output_format {
            OutputFormat::AbundanceMatrix {
                normalized,
                breakpoints,
            } => Self::AbundanceMatrix {
                normalized: Self::add_kmer_count(normalized, bf_dir),
                breakpoints,
            },
            OutputFormat::Colored { normalized } => Self::Colored {
                normalized: Self::add_kmer_count(normalized, bf_dir),
            },
            OutputFormat::Median { normalized } => Self::Median {
                normalized: Self::add_kmer_count(normalized, bf_dir),
            },
        }
    }
}
