use itertools::Itertools;
use pelt_reindeer2::pelt;

use super::super::KmerCountsAndNormalizeValue;
use super::super::compute_median;
use super::super::count_to_string_with_star;
use crate::reindeer2::query::ApproxAbundance;
use crate::reindeer2::query::format::count_to_string_with_star_normalized;

/// Computes a single cell of the matrix when we want the average of elements.
/// Returns nan if `abund_values` is empty or if it has no valid values.
pub fn cell_compute_average(abund_values: &[ApproxAbundance]) -> String {
    if abund_values.is_empty() {
        return f64::NAN.to_string();
    }
    let (sum, count) = abund_values
        .iter()
        .filter_map(|abund_values| abund_values.to_value())
        .fold((0usize, 0usize), |(s, c), abund_value| {
            (s + abund_value as usize, c + 1)
        });
    if count == 0 {
        return "/".to_string();
    }
    let average: f64 = sum as f64 / count as f64;
    average.to_string()
}

/// Computes a single cell of the matrix when we want the *normalized* average of elements.
pub fn cell_compute_average_normalized(
    abund_values: &[ApproxAbundance],
    normalize: &KmerCountsAndNormalizeValue,
    color_id: usize,
) -> String {
    if abund_values.is_empty() {
        return f64::NAN.to_string();
    }
    let (sum, count) = abund_values
        .iter()
        .filter_map(|abund_values| abund_values.to_value())
        .fold((0usize, 0usize), |(s, c), abund_value| {
            (s + abund_value as usize, c + 1)
        });
    if count == 0 {
        return "/".to_string();
    }
    let average: f64 = sum as f64 / count as f64;
    // TODO discuss if normalizing before or after average
    let KmerCountsAndNormalizeValue {
        kmer_counts,
        normalize_value,
    } = normalize;
    (average / kmer_counts[color_id] as f64 * (*normalize_value as f64)).to_string()
}

/// Computes a single cell of the matrix when we want the median of elements.
pub fn cell_compute_median(abund_values: &[ApproxAbundance]) -> String {
    // TODO discuss do we cast in f64, just as when normalizing ?
    let valid_values = abund_values
        .iter()
        .filter_map(|value| value.to_value())
        .collect_vec();
    if !abund_values.is_empty() && valid_values.is_empty() {
        // we used to have values, but they are all gone due to the filtration :(
        return "/".to_string();
    }
    compute_median(&valid_values).to_string()
}

/// Computes a single cell of the matrix when we want the *normalized* median of elements.
pub fn cell_compute_median_normalized(
    abund_values: &[ApproxAbundance],
    normalize: &KmerCountsAndNormalizeValue,
    color_id: usize,
) -> String {
    let valid_values = abund_values
        .iter()
        .filter_map(|value| value.to_value())
        .collect_vec();
    let median = compute_median(&valid_values);
    if !abund_values.is_empty() && valid_values.is_empty() {
        // we used to have values, but they are all gone due to the filtration :(
        return "/".to_string();
    }
    // TODO discuss if normalizing before or after median
    // TODO discuss if we convert to f64 even if no normalization
    let KmerCountsAndNormalizeValue {
        kmer_counts,
        normalize_value,
    } = normalize;
    (median / kmer_counts[color_id] as f64 * (*normalize_value as f64)).to_string()
}

/// Computes a single cell of the matrix when we want an output like RD1.
pub fn cell_compute_reindeer1(abund_values: &[ApproxAbundance]) -> String {
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
    abund_values: &[ApproxAbundance],
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
pub fn cell_compute_breakpoints(abund_values: &[ApproxAbundance], penalty: f64) -> String {
    // TODO make it impossible to trigger this from the command line
    #[cfg(debug_assertions)]
    {
        let not_queried = abund_values.iter().filter(|x| !x.is_queried()).count();
        assert!(
            not_queried == 0,
            "Cannot compute breakpoints if not querying all k-mers. Please disable sampling."
        );
    }
    let abund_values = abund_values
        .iter()
        .copied()
        .map(|x| {
            u64::from(
                x.to_value()
                    .unwrap_or_else(|| panic!("could not compute the value of {:?}", x)),
            )
        }) // TODO collect all valid values before
        .collect_vec();
    let breakpoints = pelt(&abund_values, pelt_reindeer2::score, penalty);
    let breakpoints = breakpoints.iter().map(usize::to_string).join(",");
    breakpoints
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;

    use super::ApproxAbundance as Approx;

    #[test]
    fn test_cell_compute_average() {
        assert_eq!(cell_compute_average(&[Approx::new(1)]), "1");
        assert_eq!(
            cell_compute_average(&[Approx::new(1), Approx::new(2)]),
            "1.5"
        );
        assert_eq!(cell_compute_average(&[Approx::new(1), Approx::new(9)]), "5");
        assert_eq!(
            cell_compute_average(&(0..100).map(Approx::new).collect_vec()),
            "49.5",
        );
    }

    #[test]
    fn test_cell_compute_average_normalized() {
        let normalize = KmerCountsAndNormalizeValue {
            kmer_counts: vec![10, 5],
            normalize_value: 100,
        };
        assert_eq!(
            cell_compute_average_normalized(&[Approx::new(1)], &normalize, 0),
            "10",
        );
        assert_eq!(
            cell_compute_average_normalized(&[Approx::new(4)], &normalize, 1),
            "80",
        );
        assert_eq!(
            cell_compute_average_normalized(
                &[Approx::new(1), Approx::new(2), Approx::new(1)],
                &normalize,
                0,
            ),
            "13.333333333333334",
        );
    }

    #[test]
    fn test_cell_compute_median() {
        assert_eq!(cell_compute_median(&[Approx::new(1)]), "1");
        assert_eq!(
            cell_compute_median(&[Approx::new(1), Approx::new(2)]),
            "1.5",
        );
        assert_eq!(cell_compute_median(&[Approx::new(1), Approx::new(9)]), "5",);
        assert_eq!(
            cell_compute_median(&(0..100).map(Approx::new).collect_vec()),
            "49.5",
        );
    }

    #[test]
    fn test_cell_compute_median_normalized() {
        let normalize = KmerCountsAndNormalizeValue {
            kmer_counts: vec![10, 5],
            normalize_value: 100,
        };
        assert_eq!(
            cell_compute_median_normalized(&[Approx::new(1)], &normalize, 0),
            "10",
        );
        assert_eq!(
            cell_compute_median_normalized(&[Approx::new(1), Approx::new(2)], &normalize, 1),
            "30",
        );
        assert_eq!(
            cell_compute_median_normalized(&[Approx::new(1), Approx::new(9)], &normalize, 0),
            "50",
        );
        assert_eq!(
            cell_compute_median_normalized(&(0..100).map(Approx::new).collect_vec(), &normalize, 0),
            "495",
        );
    }

    #[test]
    fn test_cell_compute_reindeer1() {
        // single number
        assert_eq!(cell_compute_reindeer1(&[Approx::new(1)]), "0:1");
        // two numbers
        assert_eq!(
            cell_compute_reindeer1(&[Approx::new(1), Approx::new(2)]),
            "0:1,1:2"
        );
        assert_eq!(
            cell_compute_reindeer1(&[Approx::new(1), Approx::new(9)]),
            "0:1,1:9"
        );
        // stretch in the middle
        assert_eq!(
            cell_compute_reindeer1(&[
                Approx::new(3),
                Approx::new(4),
                Approx::new(5),
                Approx::new(1),
                Approx::new(1),
                Approx::new(1),
                Approx::new(8)
            ]),
            "0:3,1:4,2:5,3-5:1,6:8"
        );
        // stretch as start and stop
        assert_eq!(
            cell_compute_reindeer1(&[
                Approx::new(3),
                Approx::new(3),
                Approx::new(4),
                Approx::new(5),
                Approx::new(1),
                Approx::new(1),
                Approx::new(1),
                Approx::new(8),
                Approx::new(7),
                Approx::new(7),
                Approx::new(7),
                Approx::new(7)
            ]),
            "0-1:3,2:4,3:5,4-6:1,7:8,8-11:7"
        );
        // one stretch
        assert_eq!(cell_compute_reindeer1(&[Approx::new(5); 44]), "0-43:5");
    }

    #[test]
    fn test_cell_compute_reindeer1_normalized() {
        let normalize = KmerCountsAndNormalizeValue {
            kmer_counts: vec![10, 5],
            normalize_value: 100,
        };
        // single number
        assert_eq!(
            cell_compute_reindeer1_normalized(&[Approx::new(1)], &normalize, 0),
            "0:10"
        );
        // two numbers
        assert_eq!(
            cell_compute_reindeer1_normalized(&[Approx::new(1), Approx::new(2)], &normalize, 0),
            "0:10,1:20"
        );
        assert_eq!(
            cell_compute_reindeer1_normalized(&[Approx::new(1), Approx::new(9)], &normalize, 1),
            "0:20,1:180"
        );
        // stretch in the middle
        assert_eq!(
            cell_compute_reindeer1_normalized(
                &[
                    Approx::new(3),
                    Approx::new(4),
                    Approx::new(5),
                    Approx::new(1),
                    Approx::new(1),
                    Approx::new(1),
                    Approx::new(8)
                ],
                &normalize,
                0
            ),
            "0:30,1:40,2:50,3-5:10,6:80"
        );
        // stretch as start and stop
        assert_eq!(
            cell_compute_reindeer1_normalized(
                &[
                    Approx::new(3),
                    Approx::new(3),
                    Approx::new(4),
                    Approx::new(5),
                    Approx::new(1),
                    Approx::new(1),
                    Approx::new(1),
                    Approx::new(8),
                    Approx::new(7),
                    Approx::new(7),
                    Approx::new(7),
                    Approx::new(7)
                ],
                &normalize,
                0
            ),
            "0-1:30,2:40,3:50,4-6:10,7:80,8-11:70"
        );
        // one stretch
        assert_eq!(
            cell_compute_reindeer1_normalized(&[Approx::new(5); 44], &normalize, 0),
            "0-43:50"
        );
    }

    #[test]
    fn test_cell_compute_breakpoints() {
        let input = [0, 0, 0, 0, 4, 4, 4, 4, 80, 80, 80, 0, 0, 0];
        let input = input.into_iter().map(Approx::new).collect_vec();
        let penalty = 1.0;
        let raw_result = cell_compute_breakpoints(&input, penalty);
        let numbers: HashSet<usize> = raw_result
            .split(",")
            .map(|val| val.parse().unwrap())
            .collect();
        let expected = HashSet::from_iter([0, 4, 8, 11]);
        assert_eq!(numbers, expected);
    }
}
