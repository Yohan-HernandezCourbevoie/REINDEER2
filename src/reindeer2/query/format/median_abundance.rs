use bio::io::fasta;
use itertools::Itertools;
use std::io::{self, Write};

use super::KmerCountsAndNormalizeValue;
use crate::reindeer2::query::{
    ApproxAbundance,
    format::{count_to_string_witout_star_maybe_normalized, get_full_header},
};

use super::compute_median;

fn get_non_zero_value(approx_abundance: &ApproxAbundance) -> Option<u16> {
    if let Some(value) = approx_abundance.to_value()
        && value > 0
    {
        return Some(value);
    }
    None
}

fn is_zero(approx_abundance: &ApproxAbundance) -> bool {
    if let Some(value) = approx_abundance.to_value() {
        value == 0
    } else {
        false
    }
}

/// Write abundances per kmer like RD1.
pub fn write_median_abundance(
    sequence_results: &[Vec<Vec<ApproxAbundance>>],
    batch: &[fasta::Record],
    normalize: &Option<KmerCountsAndNormalizeValue>,
    coverage: f32,
    filenames: &[String],
    writer: &mut impl Write,
) -> io::Result<()> {
    // Compute medians for each sequence and each color, then write them out
    for (color_vectors, record) in sequence_results.iter().zip(batch) {
        let full_header = get_full_header(record);
        for (color_idx, abund_values) in color_vectors.iter().enumerate() {
            let filename = &filenames[color_idx];
            if !abund_values.is_empty() {
                let non_zero_values = abund_values
                    .iter()
                    .filter_map(get_non_zero_value)
                    .collect_vec();
                // OPTIMIZE maybe we can repalce the computation of zero_count, if we treat the errors as 0
                let zeros_count = abund_values
                    .iter()
                    .filter(|approx_abundance| is_zero(approx_abundance))
                    .count();
                if !non_zero_values.is_empty()
                    && (((zeros_count as f32) / (abund_values.len() as f32)) < coverage)
                {
                    let median = compute_median(&non_zero_values);
                    if median > 0.0 {
                        let median = count_to_string_witout_star_maybe_normalized(
                            median, normalize, color_idx,
                        );
                        writeln!(writer, "{},{},{}", full_header, filename, median)?;
                    }
                }
            }
        }
    }
    Ok(())
}
