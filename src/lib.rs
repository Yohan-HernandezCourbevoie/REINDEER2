//! REINDEER 2 is an efficient and scalable k-mer abundance index.
//!
//! REINDEER 2 indexes unitig files, typically obtained either:
//! - A/ by running [GGCAT](https://github.com/algbio/ggcat)
//! - B/ via the [Logan project](https://github.com/IndexThePlanet/Logan), by following [these steps](https://github.com/IndexThePlanet/Logan/blob/main/Accessions.md)   

// deny the use of unwrap: use expect instead
#![cfg_attr(not(test), deny(clippy::unwrap_used))]
// #![deny(dead_code)]
// #![deny(unused)]
#![deny(clippy::allow_attributes_without_reason)]
#![warn(clippy::missing_const_for_fn)]
mod reindeer2;

#[cfg(feature = "self-destruct")]
pub use reindeer2::FailIndexation;

pub use reindeer2::{
    BreakpointsNormalize, BuildArgs, MatrixFormat, OutputFormat, Parameters, Reindeer2,
    ReplaceOutcome, SamplingStrategy, compute_base, process_fasta_in_batches, read_fof_file,
};

pub mod query {
    pub use super::reindeer2::query::format;
    pub use super::reindeer2::query::write_kmer_query;
}

pub mod values {
    pub use super::reindeer2::{
        approximate_value, compute_base, compute_log_abundance, query::ApproxAbundance,
    };
}

pub mod merge {
    pub use super::reindeer2::merge_multiple_indexes;
}
