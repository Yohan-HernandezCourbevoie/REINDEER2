use bio::io::fasta;
use std::io::{self, Write};

use super::KmerCountsAndNormalizeValue;
use crate::reindeer2::query::format::{
    count_to_string_witout_star_maybe_normalized, get_full_header,
};

use super::compute_median;

/// Write abundances per kmer like RD1.
pub fn write_median_abundance(
    sequence_results: &[Vec<Vec<u16>>],
    batch: &[fasta::Record],
    normalize: &Option<KmerCountsAndNormalizeValue>,
    coverage: f32,
    writer: &mut impl Write,
) -> io::Result<()> {
    // Compute medians for each sequence and each color, then write them out
    for (color_vectors, record) in sequence_results.iter().zip(batch) {
        let full_header = get_full_header(record);
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
                    if median > 0.0 {
                        let median = count_to_string_witout_star_maybe_normalized(
                            median, normalize, color_idx,
                        );
                        writeln!(writer, "{},{},{}", full_header, color_idx, median)?;
                    }
                }
            }
        }
    }
    Ok(())
}
