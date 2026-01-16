use itertools::Itertools;
use pelt::pelt;

use super::super::compute_median;
use super::super::count_to_string_with_star;
use super::super::KmerCountsAndNormalizeValue;
use crate::reindeer2::query::format::count_to_string_with_star_normalized;

/// Computes a single cell of the matrix when we want the average of elements.
pub fn cell_compute_average(abund_values: &[u16]) -> f64 {
    let sum: usize = abund_values.iter().copied().map(usize::from).sum();
    let average: f64 = sum as f64 / abund_values.len() as f64;
    average
}

/// Computes a single cell of the matrix when we want the *normalized* average of elements.
pub fn cell_compute_average_normalized(
    abund_values: &[u16],
    normalize: &KmerCountsAndNormalizeValue,
    color_id: usize,
) -> f64 {
    let average = cell_compute_average(abund_values);
    // TODO discuss if normalizing before or after average
    let KmerCountsAndNormalizeValue {
        kmer_counts,
        normalize_value,
    } = normalize;
    average / kmer_counts[color_id] as f64 * (*normalize_value as f64)
}

/// Computes a single cell of the matrix when we want the median of elements.
pub fn cell_compute_median(abund_values: &[u16]) -> f64 {
    // TODO discuss do we cast in f64, just as when normalizing ?
    compute_median(abund_values) as f64
}

/// Computes a single cell of the matrix when we want the *normalized* median of elements.
pub fn cell_compute_median_normalized(
    abund_values: &[u16],
    normalize: &KmerCountsAndNormalizeValue,
    color_id: usize,
) -> f64 {
    let median = compute_median(abund_values);
    // TODO discuss if normalizing before or after median
    // TODO discuss if we convert to f64 even if no normalization
    let KmerCountsAndNormalizeValue {
        kmer_counts,
        normalize_value,
    } = normalize;
    median as f64 / kmer_counts[color_id] as f64 * (*normalize_value as f64)
}

/// Computes a single cell of the matrix when we want an output like RD1.
pub fn cell_compute_reindeer1(abund_values: &[u16]) -> String {
    let mut start = 0;
    let mut current = abund_values[0];

    (1..=abund_values.len())
        .filter_map(|i| {
            if i == abund_values.len() || abund_values[i] != current {
                // new different value or end of query => we must write
                let val_str = count_to_string_with_star(current);
                let val_str = if start + 1 == i {
                    // only one value
                    format!("{}:{}", start, val_str)
                } else {
                    // multiple values
                    format!("{}-{}:{}", start, i - 1, val_str)
                };
                if i < abund_values.len() {
                    // not the end of query
                    start = i;
                    current = abund_values[i];
                }
                Some(val_str)
            } else {
                // no new value => we do nothing
                None
            }
        })
        .join(",")
}

/// Computes a single cell of the matrix when we want a *normalized* output like RD1.
pub fn cell_compute_reindeer1_normalized(
    abund_values: &[u16],
    normalize: &KmerCountsAndNormalizeValue,
    color_id: usize,
) -> String {
    let mut start = 0;
    let mut current = abund_values[0];

    (1..=abund_values.len())
        .filter_map(|i| {
            if i == abund_values.len() || abund_values[i] != current {
                // new different value or end of query => we must write
                let val_str = count_to_string_with_star_normalized(current, normalize, color_id);
                let val_str = if start + 1 == i {
                    // only one value
                    format!("{}:{}", start, val_str)
                } else {
                    // multiple values
                    format!("{}-{}:{}", start, i - 1, val_str)
                };
                if i < abund_values.len() {
                    // not the end of query
                    start = i;
                    current = abund_values[i];
                }
                Some(val_str)
            } else {
                // no new value => we do nothing
                None
            }
        })
        .join(",")
}

/// Computes a single cell of the matrix when we want to output the breakpoints along the query.
pub fn cell_compute_breakpoints(abund_values: &[u16], penalty: f64) -> String {
    let abund_values = abund_values.iter().copied().map(u64::from).collect_vec();
    let breakpoints = pelt(&abund_values, pelt::score, penalty);
    let breakpoints = breakpoints.iter().map(usize::to_string).join(",");
    breakpoints
}
