//! REINDEER 2 is an efficient and scalable k-mer abundance index.
//!
//! REINDEER 2 indexes unitig files, typically obtained either:
//! - A/ by running [GGCAT](https://github.com/algbio/ggcat)
//! - B/ via the [Logan project](https://github.com/IndexThePlanet/Logan), by following [these steps](https://github.com/IndexThePlanet/Logan/blob/main/Accessions.md)   

// deny the use of unwrap: use expect instead
#![cfg_attr(not(test), deny(clippy::unwrap_used))]
#![deny(dead_code)]
#![deny(unused)]
#![deny(clippy::allow_attributes_without_reason)]
#![warn(clippy::missing_const_for_fn)]
// #![deny(missing_docs)]

mod reindeer2;

#[cfg(feature = "self-destruct")]
pub use reindeer2::FailIndexation;

pub use reindeer2::{
    BreakpointsNormalize, BuildArgs, MatrixFormat, OutputFormat, Parameters, Reindeer2,
    ReplaceOutcome, SamplingStrategy, compute_base, process_fasta_in_batches, read_fof_file,
};

/// Module related to the query of REINDEER2 indexes.
pub mod query {
    pub use super::reindeer2::query::format;
    pub use super::reindeer2::query::write_kmer_query;
}

/// Module related to the approximation of abundance value performed by REINDEER2 indexation and query.
///
/// REINDEER2 encode abundance values in a way that guarentees precision (see the paper for more details).
/// The functions exposed in this module allows to encode and decode abundances.
pub mod values {
    pub use super::reindeer2::{compute_base, query::ApproxAbundance};
}

/// Module related to the merging of multiple REINDEER2 indexes. The indexes must have been built with the same set of parameters (except for the input file of file).
pub mod merge {
    pub use super::reindeer2::merge_multiple_indexes;
}
