mod abundance_matrix;
mod enriched_output_format;
mod graph_coloring;
mod median_abundance;

use bio::io::fasta;
use std::io::{self, Write};

use abundance_matrix::{write_abundance_matrix, write_header_matrix};
pub use enriched_output_format::{EnrichedOutputFormat, KmerCountsAndNormalizeValue};
use graph_coloring::graph_coloring;
use median_abundance::write_median_abundance;

pub use abundance_matrix::SEPARATOR as MATRIX_SEPARATOR;

use crate::reindeer2::query::ApproxAbundance;

fn count_to_string_with_star_normalized(
    count: ApproxAbundance,
    normalize: &KmerCountsAndNormalizeValue,
    color_id: usize,
) -> String {
    let KmerCountsAndNormalizeValue {
        kmer_counts,
        normalize_value,
    } = normalize;

    if let Some(value) = count.to_value() {
        if value == 0 {
            String::from("*")
        } else {
            let normalized_count =
                value as f64 / kmer_counts[color_id] as f64 * (*normalize_value as f64);
            // TODO discuss: might be zero, but only after normlization
            normalized_count.to_string()
        }
    } else {
        if !count.is_queried() {
            String::from("/")
        } else {
            String::from("/") // TODO discuss
        }
    }
}

fn count_to_string_with_star(count: ApproxAbundance) -> String {
    if let Some(value) = count.to_value() {
        if value == 0 {
            String::from("*")
        } else {
            value.to_string()
        }
    } else {
        if !count.is_queried() {
            String::from("/")
        } else {
            String::from("/") // TODO discuss
        }
    }
}

fn count_to_string_witout_star_maybe_normalized(
    count: f64,
    normalize: &Option<KmerCountsAndNormalizeValue>,
    color_id: usize,
) -> f64 {
    match normalize {
        Some(KmerCountsAndNormalizeValue {
            kmer_counts,
            normalize_value,
        }) => count / kmer_counts[color_id] as f64 * (*normalize_value as f64),
        None => count,
    }
}

fn compute_median(values: &[u16]) -> f64 {
    compute_median_vec(values.to_vec())
}

// Like `compute_median`, but more efficient if you alread have a vec
fn compute_median_vec(mut values: Vec<u16>) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    if values.len() == 1 {
        return values[0] as f64;
    }

    values.sort_unstable();

    if values.iter().all(|&x| x == 0) {
        0.0
    } else {
        let mid = values.len() / 2;
        if values.len() % 2 == 1 {
            values[mid] as f64
        } else {
            (values[mid - 1] as f64 + values[mid] as f64) / 2.0
        }
    }
}

fn get_full_header(record: &fasta::Record) -> String {
    let id = record.id();
    let desc = record.desc().unwrap_or("");

    format!(">{} {}", id, desc).trim().to_string()
}
pub fn write_header(
    index_file_names: &[String],
    output_format: &EnrichedOutputFormat,
    mut writer: &mut impl Write,
) -> std::io::Result<()> {
    match output_format {
        EnrichedOutputFormat::Median { normalized: _ } => {
            writeln!(writer, "header,file,abundance")?;
        }
        EnrichedOutputFormat::AbundanceMatrix { format: _ } => {
            write_header_matrix(&mut writer, index_file_names, MATRIX_SEPARATOR)?;
        }
        EnrichedOutputFormat::Colored { normalized: _ } => {
            // no header
        }
    }
    Ok(())
}

pub fn write_kmer_query(
    batch: &[fasta::Record],
    output_format: &EnrichedOutputFormat,
    coverage: f32,
    sequence_results: &[Vec<Vec<ApproxAbundance>>],
    filenames: &[String],
    mut writer: &mut impl Write,
) -> io::Result<()> {
    match output_format {
        EnrichedOutputFormat::Colored { normalized } => {
            // TODO lrobidou discuss what this comment means:
            // If your graph coloring wants to read from `sequence_results`:
            // Flush the writer to separate batch outputs if needed
            graph_coloring(sequence_results, batch, normalized, writer)
        }
        EnrichedOutputFormat::AbundanceMatrix { format } => {
            write_abundance_matrix(sequence_results, batch, format, &mut writer)
        }
        EnrichedOutputFormat::Median { normalized } => write_median_abundance(
            sequence_results,
            batch,
            normalized,
            coverage,
            filenames,
            &mut writer,
        ),
    }?;
    writer.flush()?;
    Ok(())
}
