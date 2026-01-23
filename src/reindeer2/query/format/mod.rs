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

fn count_to_string_with_star_normalized(
    count: u16,
    normalize: &KmerCountsAndNormalizeValue,
    color_id: usize,
) -> String {
    let KmerCountsAndNormalizeValue {
        kmer_counts,
        normalize_value,
    } = normalize;
    let normalized_count = count as f64 / kmer_counts[color_id] as f64 * (*normalize_value as f64);
    if normalized_count == 0.0 {
        String::from("*")
    } else {
        normalized_count.to_string()
    }
}

fn count_to_string_with_star(count: u16) -> String {
    if count == 0 {
        // not normalized and 0
        String::from("*")
    } else {
        // not normalized and not 0
        count.to_string()
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
    let mut abund_sorted = values.to_vec();
    abund_sorted.sort_unstable();
    if abund_sorted.iter().all(|&x| x == 0) {
        0.0
    } else if abund_sorted.len() == 1 {
        abund_sorted[0] as f64
    } else {
        let mid = abund_sorted.len() / 2;
        if abund_sorted.len() % 2 == 1 {
            abund_sorted[mid] as f64
        } else {
            (abund_sorted[mid - 1] as f64 + abund_sorted[mid] as f64) / 2.0
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
    sequence_results: &[Vec<Vec<u16>>],
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
