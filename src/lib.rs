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
mod reindeer2;

pub use reindeer2::{
    BreakpointsNormalize, MatrixFormat, OutputFormat, Parameters, Reindeer2, ReplaceOutcome,
    SamplingStrategy, read_fof_file,
};

pub mod query {
    pub use super::reindeer2::query;
}

pub mod merge {
    pub use super::reindeer2::merge_multiple_indexes;
}
