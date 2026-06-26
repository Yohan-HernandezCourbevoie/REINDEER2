use std::io::Write;
use std::sync::atomic::Ordering;
use std::{path::Path, sync::atomic::AtomicU64};

pub fn load_atomics_from_disk(
    save_path: &Path,
    last_chunk_done: Option<usize>,
) -> (
    std::sync::atomic::AtomicU64,
    std::sync::atomic::AtomicU64,
    std::sync::atomic::AtomicU64,
) {
    let last_chunk_done = match last_chunk_done {
        None => return (AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)),
        Some(x) => x,
    };

    let total_kmers =
        std::fs::read_to_string(save_path.join(format!("total_kmers_chunk{last_chunk_done}")))
            .expect("could not read debug atomics from the save file");
    let atomic_dense_kmers_count = std::fs::read_to_string(
        save_path.join(format!("atomic_dense_kmers_count_chunk{last_chunk_done}")),
    )
    .expect("could not read debug atomics from the save file");
    let atomic_sparse_kmers_count = std::fs::read_to_string(
        save_path.join(format!("atomic_sparse_kmers_count_chunk{last_chunk_done}")),
    )
    .expect("could not read debug atomics from the save file");

    let total_kmers: u64 = total_kmers.parse().expect("");
    let atomic_dense_kmers_count: u64 = atomic_dense_kmers_count.parse().expect("");
    let atomic_sparse_kmers_count: u64 = atomic_sparse_kmers_count.parse().expect("");

    (
        AtomicU64::new(total_kmers),
        AtomicU64::new(atomic_dense_kmers_count),
        AtomicU64::new(atomic_sparse_kmers_count),
    )
}

pub fn store_atomics_to_disk(
    save_path: &Path,
    last_chunk_done: usize,
    atomic_total_kmers: &AtomicU64,
    atomic_dense_kmers_count: &AtomicU64,
    atomic_sparse_kmers_count: &AtomicU64,
) {
    let total_kmers = atomic_total_kmers.load(Ordering::SeqCst);
    let dense_kmers_count = atomic_dense_kmers_count.load(Ordering::SeqCst);
    let sparse_kmers_count = atomic_sparse_kmers_count.load(Ordering::SeqCst);

    let total_kmers = total_kmers.to_string();
    let dense_kmers_count = dense_kmers_count.to_string();
    let sparse_kmers_count = sparse_kmers_count.to_string();

    let mut f =
        std::fs::File::create(save_path.join(format!("atomic_total_kmers_chunk{last_chunk_done}")))
            .expect("Should be able to create file");
    write!(f, "{}", total_kmers)
        .expect("should have been able to write the total number of kmer for this chunk");

    let mut f = std::fs::File::create(
        save_path.join(format!("atomic_dense_kmers_count_chunk{last_chunk_done}")),
    )
    .expect("Should be able to create file");
    write!(f, "{}", dense_kmers_count)
        .expect("should have been able to write the  number of dense kmer for this chunk");

    let mut f = std::fs::File::create(
        save_path.join(format!("atomic_sparse_kmers_count_chunk{last_chunk_done}")),
    )
    .expect("Should be able to create file");
    write!(f, "{}", sparse_kmers_count)
        .expect("should have been able to write the  number of sparse kmer for this chunk");
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

        let atomic_total_kmers = AtomicU64::new(85);
        let atomic_dense_kmers_count = AtomicU64::new(89485);
        let atomic_sparse_kmers_count = AtomicU64::new(7523);

        store_atomics_to_disk(
            save_dir,
            8,
            &atomic_total_kmers,
            &atomic_dense_kmers_count,
            &atomic_sparse_kmers_count,
        );
        let (a, b, c) = load_atomics_from_disk(save_dir, Some(8));

        assert_eq!(
            atomic_total_kmers.load(Ordering::SeqCst),
            a.load(Ordering::SeqCst)
        );
        assert_eq!(
            atomic_dense_kmers_count.load(Ordering::SeqCst),
            b.load(Ordering::SeqCst)
        );

        assert_eq!(
            atomic_sparse_kmers_count.load(Ordering::SeqCst),
            c.load(Ordering::SeqCst)
        );
    }
}
