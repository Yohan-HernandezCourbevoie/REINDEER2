mod abundance_matrix;
mod graph_coloring;
mod median_abundance;

use bio::io::fasta;
use std::collections::HashMap;
use std::io::{self, Write};

use crate::reindeer2::OutputFormat;
use abundance_matrix::write_abundance_matrix;
use graph_coloring::graph_coloring;
use median_abundance::write_median_abundance;

fn count_to_string(count: u16, normalize: bool, kmer_counts: &[usize], color_id: usize) -> String {
    if normalize {
        let normalized_count = count as f64 / kmer_counts[color_id] as f64 * 1_000_000f64;
        if normalized_count == 0.0 {
            String::from("*")
        } else {
            normalized_count.to_string()
        }
    } else if count == 0 {
        // not normalized and 0
        String::from("*")
    } else {
        // not normalized and not 0
        count.to_string()
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

pub fn write_kmer_query(
    bf_dir: &str,
    color_number: usize,
    batch: &[fasta::Record],
    output_format: OutputFormat,
    coverage: f32,
    sequence_results: HashMap<usize, Vec<Vec<u16>>>,
    mut writer: &mut impl Write,
) -> io::Result<()> {
    match output_format {
        OutputFormat::Colored { normalized } => {
            // TODO lrobidou discuss what this comment means:
            // If your graph coloring wants to read from `sequence_results`:
            // Flush the writer to separate batch outputs if needed
            graph_coloring(
                bf_dir,
                &sequence_results,
                batch,
                color_number,
                normalized,
                writer,
            )
        }
        OutputFormat::AbundanceMatrix {
            normalized,
            breakpoints,
        } => write_abundance_matrix(
            bf_dir,
            &sequence_results,
            batch,
            color_number,
            breakpoints,
            normalized,
            &mut writer,
        ),
        OutputFormat::Median { normalized } => write_median_abundance(
            bf_dir,
            &sequence_results,
            batch,
            color_number,
            normalized,
            coverage,
            &mut writer,
        ),
    }?;
    writer.flush()?;
    Ok(())
}
