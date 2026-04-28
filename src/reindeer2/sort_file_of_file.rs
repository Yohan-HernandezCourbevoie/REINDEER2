use std::fs;

use itertools::Itertools;
use rayon::prelude::*;

/// Sorts a list of file paths by the size of the files they point to (ascending).
/// Paths that cannot be stat'd are placed at the end in their original relative order.
pub fn sort_paths_by_file_size(paths: &mut Vec<String>) {
    // Gather sizes upfront to avoid redundant stat calls during sort comparisons
    let sizes: Vec<Option<u64>> = paths
        .par_iter()
        .map(|p| {
            fs::metadata(p)
                .map(|m| m.len())
                .map_err(|e| eprintln!("warn: cannot stat {p:?}: {e}"))
                .ok()
        })
        .collect();

    let mut indexed = sizes.into_iter().enumerate().collect_vec();

    indexed.sort_by(|(_, a), (_, b)| match (a, b) {
        (Some(sa), Some(sb)) => sa.cmp(sb),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => std::cmp::Ordering::Equal,
    });

    let original = std::mem::take(paths);
    *paths = indexed
        .into_iter()
        .map(|(i, _)| original[i].clone())
        .collect();
}
