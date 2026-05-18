//! This module handles writing the matrix to a writer.
//! We support multiple formats. To prevent code duplication, two functions
//! (namely `write_matrix_content` and `write_normalized_matrix_content`)
//! are provided.
//!
//! These functions are used as a "templates" to generate specialized matrix writing functions.
//! Each "template" accepts a function that computes a single cell of the matrix
//! (these functions are in the module `cell_computation`),
//! and uses it to compute the whole matrix.
//! A nice bonus of this design is that
//! the test to decide how to output the matrix is done only once, when instantiating the template
//! (as opposed to in the middle of the hot loop),
//! opening the door to optimization from the compiler.

mod cell_computation;

pub const SEPARATOR: &str = "\t";

use std::fmt::Display;
use std::io::{self, Write};

use bio::io::fasta;

use super::KmerCountsAndNormalizeValue;
use crate::reindeer2::query::ApproxAbundance;
use crate::reindeer2::query::format::enriched_output_format::{
    BreakpointsXorEnrichedNormalize, EnrichedMatrixFormat,
};

use cell_computation::{
    cell_compute_average, cell_compute_average_normalized, cell_compute_breakpoints,
    cell_compute_median, cell_compute_median_normalized, cell_compute_reindeer1,
    cell_compute_reindeer1_normalized,
};

/// Writes the header of the matrix, including a `\n` at the end
pub fn write_header_matrix(
    writer: &mut impl Write,
    indexed_files: &[String],
    sep: &str,
) -> io::Result<()> {
    write!(writer, "query")?;
    for indexed_file in indexed_files {
        write!(writer, "{sep}{indexed_file}")?;
    }
    writeln!(writer)?;
    Ok(())
}

/// Write the matrix to the writer according to some function.
fn write_matrix_content<F, W, R>(
    sequence_results: &[Vec<Vec<ApproxAbundance>>],
    batch: &[fasta::Record],
    cell_computation: F,
    writer: &mut W,
) -> io::Result<()>
where
    F: Fn(&[ApproxAbundance]) -> R,
    W: Write,
    R: Display,
{
    for (color_vectors, record) in sequence_results.iter().zip(batch) {
        write!(writer, "{}", record.id())?;
        for abund_values in color_vectors {
            let value = cell_computation(abund_values);
            write!(writer, "{SEPARATOR}{value}")?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

/// Write the *normalized* matrix to the writer according to some function.
fn write_normalized_matrix_content<F, W, R>(
    sequence_results: &[Vec<Vec<ApproxAbundance>>],
    batch: &[fasta::Record],
    normalize: &KmerCountsAndNormalizeValue,
    cell_computation: F,
    writer: &mut W,
) -> io::Result<()>
where
    F: Fn(&[ApproxAbundance], &KmerCountsAndNormalizeValue, usize) -> R,
    W: Write,
    R: Display,
{
    for (color_vectors, record) in sequence_results.iter().zip(batch) {
        write!(writer, "{}", record.id())?;
        for (color_id, abund_values) in color_vectors.iter().enumerate() {
            let value = cell_computation(abund_values, normalize, color_id);
            write!(writer, "{SEPARATOR}{value}")?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

// Write abundance matrix to the writer.
pub fn write_abundance_matrix(
    sequence_results: &[Vec<Vec<ApproxAbundance>>],
    batch: &[fasta::Record],
    format: &EnrichedMatrixFormat,
    writer: &mut impl Write,
) -> io::Result<()> {
    use EnrichedMatrixFormat as MatrixFormat;
    match format {
        MatrixFormat::Average { normalized } => match normalized {
            Some(normalize) => write_normalized_matrix_content(
                sequence_results,
                batch,
                normalize,
                cell_compute_average_normalized,
                writer,
            ),
            None => write_matrix_content(sequence_results, batch, cell_compute_average, writer),
        },
        MatrixFormat::Raw(raw_infos) => {
            match raw_infos {
                Some(BreakpointsXorEnrichedNormalize::Breakpoints(penalty)) => {
                    // Here, we would like to use `cell_compute_breakpoints` as an argument of `write_matrix_content`.
                    // Unfortunately, its signature doesn't match what `write_matrix_content` accepts.
                    // Therefore, we use a lambda that captures `penalty` and ignore the parameters we are not interested in.
                    let cell_compute = |abund_values: &[ApproxAbundance]| {
                        cell_compute_breakpoints(abund_values, *penalty)
                    };
                    // we can now pass the lambda to `write_matrix_content`
                    write_matrix_content(sequence_results, batch, cell_compute, writer)
                }
                Some(BreakpointsXorEnrichedNormalize::Normalize(normalized)) => {
                    write_normalized_matrix_content(
                        sequence_results,
                        batch,
                        normalized,
                        cell_compute_reindeer1_normalized,
                        writer,
                    )
                }
                None => {
                    write_matrix_content(sequence_results, batch, cell_compute_reindeer1, writer)
                }
            }
        }
        MatrixFormat::Median { normalized } => match normalized {
            Some(normalize) => write_normalized_matrix_content(
                sequence_results,
                batch,
                normalize,
                cell_compute_median_normalized,
                writer,
            ),
            None => write_matrix_content(sequence_results, batch, cell_compute_median, writer),
        },
    }
}
