mod abundance_matrix;
mod enriched_output_format;
mod graph_coloring;
mod median_abundance;

use bio::io::fasta;
use std::collections::HashMap;
use std::io::{self, Write};
use std::path::Path;

use abundance_matrix::{write_abundance_matrix, write_header_matrix};
pub use enriched_output_format::{EnrichedOutputFormat, KmerCountsAndNormalizeValue};
use graph_coloring::graph_coloring;
use median_abundance::write_median_abundance;

pub use abundance_matrix::SEPARATOR as MATRIX_SEPARATOR;

use super::read_indexed_file_names;

fn count_to_string_with_star(
    count: u16,
    normalize: &Option<KmerCountsAndNormalizeValue>,
    color_id: usize,
) -> String {
    match normalize {
        Some(KmerCountsAndNormalizeValue {
            kmer_counts,
            normalize_value,
        }) => {
            let normalized_count =
                count as f64 / kmer_counts[color_id] as f64 * (*normalize_value as f64);
            if normalized_count == 0.0 {
                String::from("*")
            } else {
                normalized_count.to_string()
            }
        }
        None => {
            if count == 0 {
                // not normalized and 0
                String::from("*")
            } else {
                // not normalized and not 0
                count.to_string()
            }
        }
    }
}

fn count_to_string_witout_star(
    count: u16,
    normalize: &Option<KmerCountsAndNormalizeValue>,
    color_id: usize,
) -> f64 {
    match normalize {
        Some(KmerCountsAndNormalizeValue {
            kmer_counts,
            normalize_value,
        }) => count as f64 / kmer_counts[color_id] as f64 * (*normalize_value as f64),
        None => count as f64,
    }
}

fn compute_median(values: &[u16]) -> u16 {
    let mut abund_sorted = values.to_vec();
    abund_sorted.sort_unstable();
    if abund_sorted.iter().all(|&x| x == 0) {
        0
    } else if abund_sorted.len() == 1 {
        abund_sorted[0]
    } else {
        let mid = abund_sorted.len() / 2;
        if abund_sorted.len() % 2 == 1 {
            abund_sorted[mid]
        } else {
            (abund_sorted[mid - 1] + abund_sorted[mid]) / 2
        }
    }
}

fn get_full_header(record: &fasta::Record) -> String {
    let id = record.id();
    let desc = record.desc().unwrap_or("");

    format!(">{} {}", id, desc).trim().to_string()
}
pub fn write_header(
    bf_dir: &str,
    output_format: &EnrichedOutputFormat,
    mut writer: &mut impl Write,
) -> std::io::Result<()> {
    match output_format {
        EnrichedOutputFormat::Median { normalized: _ } => {
            writeln!(writer, "header,file,abundance")?;
        }
        EnrichedOutputFormat::AbundanceMatrix {
            normalized: _,
            breakpoints: _,
        } => {
            let indexed_files = Path::new(&bf_dir).join("path.txt");
            let indexed_files: Vec<String> = read_indexed_file_names(indexed_files);
            write_header_matrix(&mut writer, indexed_files, MATRIX_SEPARATOR)?;
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
    sequence_results: HashMap<usize, Vec<Vec<u16>>>,
    mut writer: &mut impl Write,
) -> io::Result<()> {
    match output_format {
        EnrichedOutputFormat::Colored { normalized } => {
            // TODO lrobidou discuss what this comment means:
            // If your graph coloring wants to read from `sequence_results`:
            // Flush the writer to separate batch outputs if needed
            graph_coloring(&sequence_results, batch, normalized, writer)
        }
        EnrichedOutputFormat::AbundanceMatrix {
            normalized,
            breakpoints,
        } => write_abundance_matrix(
            &sequence_results,
            batch,
            breakpoints,
            normalized,
            &mut writer,
        ),
        EnrichedOutputFormat::Median { normalized } => {
            write_median_abundance(&sequence_results, batch, normalized, coverage, &mut writer)
        }
    }?;
    writer.flush()?;
    Ok(())
}
