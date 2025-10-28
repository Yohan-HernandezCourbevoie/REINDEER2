use bio::io::fasta;
use std::collections::HashMap;
use std::io::{self, Write};

use super::super::load_kmer_counts_vector;
use super::compute_median;

/// Write abundances per kmer like RD1.
pub fn write_median_abundance(
    bf_dir: &str,
    sequence_results: &HashMap<usize, Vec<Vec<u16>>>,
    batch: &[fasta::Record],
    color_number: usize,
    normalize: bool,
    coverage: f32,
    writer: &mut impl Write,
) -> io::Result<()> {
    writeln!(writer, "header,file,abundance")?;
    // we need the count of kmers if we want to normalize them
    let kmer_counts = if normalize {
        load_kmer_counts_vector(bf_dir).expect("Failed to load from disk the kmer counts vector")
    } else {
        vec![color_number, 0] // TODO bizarre
    };
    // Compute medians for each sequence and each color, then write them out
    for (record_id, color_vectors) in sequence_results {
        let record = &batch[*record_id];
        let id = record.id();
        let desc = record.desc().unwrap_or("");
        let full_header = format!(">{} {}", id, desc).trim().to_string();
        for (color_idx, abund_values) in color_vectors.iter().enumerate() {
            if !abund_values.is_empty() {
                let mut zeros_count = 0;
                let mut non_zero_values: Vec<u16> = Vec::new();
                abund_values.iter().for_each(|value| {
                    if *value == 0 {
                        zeros_count += 1;
                    } else {
                        non_zero_values.push(*value);
                    }
                });
                if !non_zero_values.is_empty()
                    && (((zeros_count as f32) / (abund_values.len() as f32)) < coverage)
                {
                    let median = compute_median(&non_zero_values);
                    if median > 0 {
                        let median = if normalize {
                            median as f64 / kmer_counts[color_idx] as f64 * 1_000_000f64
                        } else {
                            median as f64
                        };
                        writeln!(writer, "{},{},{}", full_header, color_idx, median)?;
                    }
                }
            }
        }
    }
    Ok(())
}
