#![cfg(any(debug_assertions, test))]

use std::io::Write;
use std::sync::atomic::Ordering;
use std::{path::Path, sync::atomic::AtomicU64};

pub fn load_debug_atomics_from_disk(
    save_path: &Path,
    last_chunk_done: Option<usize>,
) -> (std::sync::atomic::AtomicU64, std::sync::atomic::AtomicU64) {
    let last_chunk_done = match last_chunk_done {
        None => return (AtomicU64::new(0), AtomicU64::new(0)),
        Some(x) => x,
    };
    let atomic_sparse_one_seen = std::fs::read_to_string(
        save_path.join(format!("atomic_sparse_one_seen_chunk{last_chunk_done}")),
    )
    .expect("could not read debug atomics from the save file");
    let atomic_sparse_fp_seen = std::fs::read_to_string(
        save_path.join(format!("atomic_sparse_fp_seen_chunk{last_chunk_done}")),
    )
    .expect("could not read debug atomics from the save file");
    let atomic_sparse_one_seen: u64 = atomic_sparse_one_seen.parse().expect("");
    let atomic_sparse_fp_seen: u64 = atomic_sparse_fp_seen.parse().expect("");

    (
        AtomicU64::new(atomic_sparse_one_seen),
        AtomicU64::new(atomic_sparse_fp_seen),
    )
}

pub fn store_debug_atomics_to_disk(
    save_path: &Path,
    last_chunk_done: usize,
    atomic_sparse_one_seen: &AtomicU64,
    atomic_sparse_fp_seen: &AtomicU64,
) {
    let sparse_one_seen = atomic_sparse_one_seen.load(Ordering::SeqCst);
    let sparse_fp_seenlet = atomic_sparse_fp_seen.load(Ordering::SeqCst);

    let sparse_one_seen = sparse_one_seen.to_string();
    let sparse_fp_seen = sparse_fp_seenlet.to_string();

    let mut f = std::fs::File::create(
        save_path.join(format!("atomic_sparse_one_seen_chunk{last_chunk_done}")),
    )
    .expect("Should be able to create file");
    write!(f, "{}", sparse_one_seen)
        .expect("should have been able to write the number of ones seen for this chunk");

    let mut f = std::fs::File::create(
        save_path.join(format!("atomic_sparse_fp_seen_chunk{last_chunk_done}")),
    )
    .expect("Should be able to create file");
    write!(f, "{}", sparse_fp_seen)
        .expect("should have been able to write the number of fp seen for this chunk");
}

#[cfg(test)]
mod tests {
    use crate::reindeer2::test_utils::AutoRemoveDirectory;

    use super::*;

    use rstest::{fixture, rstest};

    #[fixture]
    pub fn random_directory() -> AutoRemoveDirectory {
        AutoRemoveDirectory::create_random()
    }

    #[rstest]
    fn test_saves(random_directory: AutoRemoveDirectory) {
        let save_dir = random_directory.filename().to_str().unwrap();
        std::fs::create_dir_all(save_dir).expect("Failed to create test directory");
        let save_dir = Path::new(save_dir);

        let atomic_sparse_one_seen = std::sync::atomic::AtomicU64::new(85);
        let atomic_sparse_fp_seen = std::sync::atomic::AtomicU64::new(89485);

        store_debug_atomics_to_disk(save_dir, 8, &atomic_sparse_one_seen, &atomic_sparse_fp_seen);
        let (a, b) = load_debug_atomics_from_disk(save_dir, Some(8));

        assert_eq!(
            atomic_sparse_one_seen.load(Ordering::SeqCst),
            a.load(Ordering::SeqCst)
        );
        assert_eq!(
            atomic_sparse_fp_seen.load(Ordering::SeqCst),
            b.load(Ordering::SeqCst)
        );
    }
}
