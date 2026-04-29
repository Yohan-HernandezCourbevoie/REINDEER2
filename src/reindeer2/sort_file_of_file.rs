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

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::tempdir;

    // Helper: create a temp file with a specific size (in bytes)
    fn make_temp_file(dir: &tempfile::TempDir, name: &str, size: u64) -> String {
        let path = dir.path().join(name);
        let data = vec![0u8; size as usize];
        fs::write(&path, data).unwrap();
        path.to_str().unwrap().to_owned()
    }

    #[test]
    fn test_sorts_ascending_by_size() {
        let dir = tempdir().unwrap();
        let small = make_temp_file(&dir, "small.bin", 100);
        let medium = make_temp_file(&dir, "medium.bin", 200);
        let large = make_temp_file(&dir, "large.bin", 300);

        let mut paths = vec![large.clone(), small.clone(), medium.clone()];
        sort_paths_by_file_size(&mut paths);

        assert_eq!(paths, vec![small, medium, large]);
    }

    #[test]
    fn test_already_sorted_unchanged() {
        let dir = tempdir().unwrap();
        let a = make_temp_file(&dir, "a.bin", 10);
        let b = make_temp_file(&dir, "b.bin", 20);
        let c = make_temp_file(&dir, "c.bin", 30);

        let mut paths = vec![a.clone(), b.clone(), c.clone()];
        sort_paths_by_file_size(&mut paths);

        assert_eq!(paths, vec![a, b, c]);
    }

    #[test]
    fn test_reverse_sorted() {
        let dir = tempdir().unwrap();
        let a = make_temp_file(&dir, "a.bin", 10);
        let b = make_temp_file(&dir, "b.bin", 20);
        let c = make_temp_file(&dir, "c.bin", 30);

        let mut paths = vec![c.clone(), b.clone(), a.clone()];
        sort_paths_by_file_size(&mut paths);

        assert_eq!(paths, vec![a, b, c]);
    }

    #[test]
    fn test_missing_paths_sorted_to_end() {
        let dir = tempdir().unwrap();
        let small = make_temp_file(&dir, "small.bin", 50);
        let large = make_temp_file(&dir, "large.bin", 150);
        let ghost = "/nonexistent/path/ghost.bin".to_owned();

        let mut paths = vec![ghost.clone(), large.clone(), small.clone()];
        sort_paths_by_file_size(&mut paths);

        assert_eq!(paths[0], small);
        assert_eq!(paths[1], large);
        assert_eq!(paths[2], ghost); // missing file pushed to the end
    }

    #[test]
    fn test_multiple_missing_paths_preserve_relative_order() {
        let dir = tempdir().unwrap();
        let real = make_temp_file(&dir, "real.bin", 100);
        let ghost1 = "/nonexistent/ghost1.bin".to_owned();
        let ghost2 = "/nonexistent/ghost2.bin".to_owned();

        let mut paths = vec![ghost1.clone(), real.clone(), ghost2.clone()];
        sort_paths_by_file_size(&mut paths);

        assert_eq!(paths[0], real);
        // ghost1 came before ghost2 in the input, order should be preserved
        assert_eq!(paths[1], ghost1);
        assert_eq!(paths[2], ghost2);
    }

    #[test]
    fn test_files_with_equal_sizes_preserve_relative_order() {
        let dir = tempdir().unwrap();
        let a = make_temp_file(&dir, "a.bin", 100);
        let b = make_temp_file(&dir, "b.bin", 100);
        let c = make_temp_file(&dir, "c.bin", 100);

        let mut paths = vec![a.clone(), b.clone(), c.clone()];
        sort_paths_by_file_size(&mut paths);

        // All equal: original order should be preserved (stable sort)
        assert_eq!(paths, vec![a, b, c]);
    }

    #[test]
    fn test_empty_input() {
        let mut paths: Vec<String> = vec![];
        sort_paths_by_file_size(&mut paths);
        assert!(paths.is_empty());
    }

    #[test]
    fn test_single_entry() {
        let dir = tempdir().unwrap();
        let a = make_temp_file(&dir, "a.bin", 42);

        let mut paths = vec![a.clone()];
        sort_paths_by_file_size(&mut paths);

        assert_eq!(paths, vec![a]);
    }

    #[test]
    fn test_empty_file() {
        let dir = tempdir().unwrap();
        let empty = make_temp_file(&dir, "empty.bin", 0);
        let small = make_temp_file(&dir, "small.bin", 10);

        let mut paths = vec![small.clone(), empty.clone()];
        sort_paths_by_file_size(&mut paths);

        assert_eq!(paths, vec![empty, small]);
    }
}
