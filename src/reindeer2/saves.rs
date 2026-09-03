//! Saves states for REINDEER2

use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::sync::{Arc, Mutex};
use std::{
    collections::HashSet,
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
};

use crate::reindeer2::NB_FILE_IN_AN_INDEX;

/// Represent the last chunk completed before a crash happened during the chunk indexation phase.
#[derive(Debug, PartialEq, Eq)]
pub struct Chunks {
    /// The id of the last chunk sucessfully completed. Can be None if the crash happened before the first chunk was completed.
    done: Option<usize>,
}

/// Represent the steps completed before a crash happened during the merge indexation phase.
#[derive(Debug, PartialEq, Eq)]
pub struct Merge {
    /// The set of merge sucessfully done before a crash happened.
    merges_done: HashSet<usize>,
}

/// Saves of the index
#[derive(Debug, PartialEq, Eq)]
pub struct Saves<T> {
    path: PathBuf,
    inner: T,
}

#[derive(Debug, PartialEq)]
pub enum SavesState {
    Chunk { saves: Saves<Chunks> },
    Merge { saves: Saves<Merge> },
}

impl<T> Saves<T> {
    pub fn get_filename_for_chunk(chunk: usize) -> String {
        format!("chunk_{chunk}_done")
    }

    pub const fn get_filename_for_all_chunk_done() -> &'static str {
        "chunks_all_done"
    }

    pub fn get_filename_for_merge(merge: usize) -> String {
        format!("merge_{merge}_done")
    }

    pub fn detect_crash(path: &Path, nb_chunk: usize) -> SavesState {
        let crashed_after_chunk =
            std::fs::exists(path.join(Self::get_filename_for_all_chunk_done()))
                .expect("should have been able to read the save folder");
        if crashed_after_chunk {
            let mut merges_done = HashSet::new();
            for i in 0..NB_FILE_IN_AN_INDEX {
                let chunk_filename = Self::get_filename_for_merge(i);
                let chunk_path = path.join(&chunk_filename);
                let chunk_done = std::fs::exists(chunk_path)
                    .expect("should have been able to check if the chunk was done");
                if chunk_done {
                    merges_done.insert(i);
                }
            }

            // TODO check we cannot fire this one
            assert!(merges_done.len() != NB_FILE_IN_AN_INDEX);

            SavesState::Merge {
                saves: Saves {
                    path: path.to_path_buf(),
                    inner: Merge { merges_done },
                },
            }
        } else {
            // we crached before finishing all the chunks
            let mut last_done = None;
            for i in 0..nb_chunk {
                let chunk_save_name = Self::get_filename_for_chunk(i);
                let chunk_done = std::fs::exists(path.join(&chunk_save_name))
                    .expect("should have been able to read the save folder");
                if chunk_done {
                    last_done = Some(i);
                }
            }

            SavesState::Chunk {
                saves: Saves {
                    path: path.to_path_buf(),
                    inner: Chunks { done: last_done },
                },
            }
        }
    }
}

impl Saves<Chunks> {
    pub const fn new(path: PathBuf) -> Self {
        Self {
            path,
            inner: Chunks { done: None },
        }
    }

    pub const fn get_last_done(&self) -> Option<usize> {
        self.inner.done
    }

    pub const fn get_current_chunk(&self) -> usize {
        match self.inner.done {
            None => 0,
            Some(x) => x + 1,
        }
    }

    pub const fn get_nb_chunk_that_can_be_skipped(&self) -> usize {
        // same as get_current_chunk
        // but I leave the code duplicated in case the two have to diverge at some point
        match self.inner.done {
            None => 0,
            Some(x) => x + 1,
        }
    }

    pub fn one_chunk_done(&mut self, chunk_done: usize) {
        match &mut self.inner.done {
            None => {
                assert!(chunk_done == 0);
                self.inner.done = Some(0);
            }
            Some(i) => {
                *i += 1;
                assert!(*i == chunk_done);
            }
        }
        let chunk_save_name = Self::get_filename_for_chunk(chunk_done);
        // TODO check the error message
        std::fs::File::create_new(self.path.join(&chunk_save_name)).unwrap_or_else(|err| {
            panic!("attempted to save a chunk that was already saved: {chunk_save_name} ({err})")
        });
        if chunk_done > 0 {
            let previous_chunk = chunk_done - 1;
            let previous_chunk_save_name = Self::get_filename_for_chunk(previous_chunk);
            std::fs::remove_file(self.path.join(&previous_chunk_save_name)).unwrap_or_else(|err| {
                panic!(
                    "failed to remove a save for the previous chunk\ncurrent chunk: {chunk_save_name}, previous_chunk: {previous_chunk_save_name} ({err})"
                )
            });
        }
    }

    pub fn chunks_all_done(self) -> Saves<Merge> {
        std::fs::File::create_new(self.path.join(Self::get_filename_for_all_chunk_done()))
            .unwrap_or_else(|err| {
                panic!("attempted to create a save after indexing all chunks but failed ({err})")
            });

        if let Some(chunk_id) = self.inner.done {
            let chunk_save_name = Self::get_filename_for_chunk(chunk_id);
            std::fs::remove_file(self.path.join(&chunk_save_name)).unwrap_or_else(|err| {
                panic!(
                    "failed to remove a save for the last chunk\nprevious_chunk: {chunk_save_name} ({err})"
                )
            })
        }

        Saves::<Merge> {
            path: self.path,
            inner: Merge {
                merges_done: HashSet::new(),
            },
        }
    }

    #[cfg(any(debug_assertions, test))]
    pub fn load_debug_atomics_from_disk(
        &self,
    ) -> (std::sync::atomic::AtomicU64, std::sync::atomic::AtomicU64) {
        let last_chunk_done = match self.get_last_done() {
            None => return (AtomicU64::new(0), AtomicU64::new(0)),
            Some(x) => x,
        };
        let atomic_sparse_one_seen = std::fs::read_to_string(
            self.path
                .join(format!("atomic_sparse_one_seen_chunk{last_chunk_done}")),
        )
        .expect("could not read debug atomics from the save file");
        let atomic_sparse_fp_seen = std::fs::read_to_string(
            self.path
                .join(format!("atomic_sparse_fp_seen_chunk{last_chunk_done}")),
        )
        .expect("could not read debug atomics from the save file");
        let atomic_sparse_one_seen: u64 = atomic_sparse_one_seen.parse().expect("");
        let atomic_sparse_fp_seen: u64 = atomic_sparse_fp_seen.parse().expect("");

        (
            AtomicU64::new(atomic_sparse_one_seen),
            AtomicU64::new(atomic_sparse_fp_seen),
        )
    }

    #[cfg(any(debug_assertions, test))]
    pub fn store_current_chunk_debug_atomics_to_disk(
        &self,
        atomic_sparse_one_seen: &AtomicU64,
        atomic_sparse_fp_seen: &AtomicU64,
    ) {
        let current_chunk = self.get_current_chunk();
        let sparse_one_seen = atomic_sparse_one_seen.load(Ordering::SeqCst);
        let sparse_fp_seenlet = atomic_sparse_fp_seen.load(Ordering::SeqCst);

        let sparse_one_seen = sparse_one_seen.to_string();
        let sparse_fp_seen = sparse_fp_seenlet.to_string();

        let mut f = std::fs::File::create(
            self.path
                .join(format!("atomic_sparse_one_seen_chunk{current_chunk}")),
        )
        .expect("Should be able to create file");
        write!(f, "{}", sparse_one_seen)
            .expect("should have been able to write the number of ones seen for this chunk");

        let mut f = std::fs::File::create(
            self.path
                .join(format!("atomic_sparse_fp_seen_chunk{current_chunk}")),
        )
        .expect("Should be able to create file");
        write!(f, "{}", sparse_fp_seen)
            .expect("should have been able to write the number of fp seen for this chunk");
    }

    pub fn load_atomics_from_disk(
        &self,
    ) -> (
        std::sync::atomic::AtomicU64,
        std::sync::atomic::AtomicU64,
        std::sync::atomic::AtomicU64,
    ) {
        let last_chunk_done = match self.get_last_done() {
            None => return (AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)),
            Some(x) => x,
        };

        let total_kmers = std::fs::read_to_string(
            self.path
                .join(format!("atomic_total_kmers_chunk{last_chunk_done}")),
        )
        .expect("could not read debug atomics from the save file");
        let atomic_dense_kmers_count = std::fs::read_to_string(
            self.path
                .join(format!("atomic_dense_kmers_count_chunk{last_chunk_done}")),
        )
        .expect("could not read debug atomics from the save file");
        let atomic_sparse_kmers_count = std::fs::read_to_string(
            self.path
                .join(format!("atomic_sparse_kmers_count_chunk{last_chunk_done}")),
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

    pub fn store_current_atomics_to_disk(
        &self,
        atomic_total_kmers: &AtomicU64,
        atomic_dense_kmers_count: &AtomicU64,
        atomic_sparse_kmers_count: &AtomicU64,
    ) {
        let current_chunk = self.get_current_chunk();
        let total_kmers = atomic_total_kmers.load(Ordering::SeqCst);
        let dense_kmers_count = atomic_dense_kmers_count.load(Ordering::SeqCst);
        let sparse_kmers_count = atomic_sparse_kmers_count.load(Ordering::SeqCst);

        let total_kmers = total_kmers.to_string();
        let dense_kmers_count = dense_kmers_count.to_string();
        let sparse_kmers_count = sparse_kmers_count.to_string();

        let mut f = std::fs::File::create(
            self.path
                .join(format!("atomic_total_kmers_chunk{current_chunk}")),
        )
        .expect("Should be able to create file");
        write!(f, "{}", total_kmers)
            .expect("should have been able to write the total number of kmer for this chunk");

        let mut f = std::fs::File::create(
            self.path
                .join(format!("atomic_dense_kmers_count_chunk{current_chunk}")),
        )
        .expect("Should be able to create file");
        write!(f, "{}", dense_kmers_count)
            .expect("should have been able to write the  number of dense kmer for this chunk");

        let mut f = std::fs::File::create(
            self.path
                .join(format!("atomic_sparse_kmers_count_chunk{current_chunk}")),
        )
        .expect("Should be able to create file");
        write!(f, "{}", sparse_kmers_count)
            .expect("should have been able to write the  number of sparse kmer for this chunk");
    }

    pub fn store_current_kmer_counts_to_disk(
        &self,
        kmer_counts_vector: &Arc<Mutex<Vec<usize>>>,
    ) -> io::Result<()> {
        let current_chunk = self.get_current_chunk();
        // OPTIMIZE we may be able to drop the lock before writing to disk
        let file_path = self.path.join(format!(
            "kmer_counts_per_color_after_chunk_{current_chunk}.bin"
        ));
        let file = File::create(&file_path)?;
        let mut writer = BufWriter::new(file);
        let locked_vector = kmer_counts_vector.lock().expect(
            "fatal error: a thread holding the mutex panicked, so this thread will panic as well",
        );
        let locked_vector_ref: &Vec<usize> = &locked_vector;
        let binary_encoded = bincode::serialize(locked_vector_ref)
            .expect("should have been able to serialize the count of k-mers");
        drop(locked_vector);
        writer.write_all(&binary_encoded)?;
        Ok(())
    }

    pub fn load_kmer_counts_vector(&self, nb_color: usize) -> io::Result<Vec<usize>> {
        let chunk = match self.get_last_done() {
            None => return Ok(vec![0; nb_color]),
            Some(chunk) => chunk,
        };

        let file_path =
            Path::new(&self.path).join(format!("kmer_counts_per_color_after_chunk_{chunk}.bin"));
        let mut file = File::open(file_path)?;

        // Read the rest of the file to deserialize the hashmap
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;
        let counts_vector = bincode::deserialize_from(&buffer[..]).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "Failed to deserialize the counts vector",
            )
        })?;
        Ok(counts_vector)
    }
}

impl Saves<Merge> {
    // panics if the merge was already done
    pub fn one_merge_done(&mut self, merge_done: usize) {
        let inserted = self.inner.merges_done.insert(merge_done);
        assert!(inserted, "attempted to merge two times the same files");

        std::fs::File::create_new(self.path.join(Self::get_filename_for_merge(merge_done)))
            .unwrap_or_else(|err| {
                panic!("attempted to create a save after merging {merge_done} but failed ({err})")
            });
    }

    pub fn merge_all_done(self) {
        for i in 0..NB_FILE_IN_AN_INDEX {
            let chunk_filename = Self::get_filename_for_merge(i);
            let chunk_path = self.path.join(&chunk_filename);
            let chunk_done = std::fs::exists(chunk_path)
                .expect("should have been able to check if the chunk was done");
            assert!(chunk_done);
        }
    }

    pub fn get_still_to_be_merged(&self) -> std::collections::HashSet<usize> {
        let mut still_to_be_merged = HashSet::new();
        for i in 0..NB_FILE_IN_AN_INDEX {
            if !self.inner.merges_done.contains(&i) {
                still_to_be_merged.insert(i);
            }
        }
        still_to_be_merged
    }

    pub fn get_merge_done(&self) -> std::collections::HashSet<usize> {
        self.inner.merges_done.clone()
    }

    pub fn from_disk_state(path: PathBuf) -> Self {
        let all_chunk_done = std::fs::exists(path.join(Self::get_filename_for_all_chunk_done()))
            .expect("should have been able to check if all chunk are done");
        assert!(
            all_chunk_done,
            "attempting to resume the merge while all chunks are not done"
        );

        let mut merges_done = HashSet::new();
        for i in 0..NB_FILE_IN_AN_INDEX {
            let chunk_filename = Self::get_filename_for_merge(i);
            let chunk_path = path.join(&chunk_filename);
            let chunk_done = std::fs::exists(chunk_path)
                .expect("should have been able to check if the chunk was done");
            if chunk_done {
                merges_done.insert(i);
            }
        }
        Self {
            path,
            inner: Merge { merges_done },
        }
    }

    pub fn nb_merge_left_to_be_done(&self) -> usize {
        NB_FILE_IN_AN_INDEX - self.inner.merges_done.len()
    }
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
        let save_dir = Path::new(save_dir);
        std::fs::create_dir_all(save_dir).expect("Failed to create test directory");

        let mut saves = Saves::new(PathBuf::from(save_dir));
        assert_eq!(saves.inner, Chunks { done: None });
        assert_eq!(
            Saves::<Chunks>::detect_crash(save_dir, 3),
            SavesState::Chunk {
                saves: Saves {
                    path: PathBuf::from(save_dir),
                    inner: Chunks { done: None }
                }
            }
        );

        saves.one_chunk_done(0);
        assert_eq!(saves.inner, Chunks { done: Some(0) });
        assert_eq!(
            Saves::<Chunks>::detect_crash(save_dir, 3),
            SavesState::Chunk {
                saves: Saves {
                    path: PathBuf::from(save_dir),
                    inner: Chunks { done: Some(0) }
                }
            }
        );

        saves.one_chunk_done(1);
        assert_eq!(saves.inner, Chunks { done: Some(1) });
        assert_eq!(
            Saves::<Chunks>::detect_crash(save_dir, 3),
            SavesState::Chunk {
                saves: Saves {
                    path: PathBuf::from(save_dir),
                    inner: Chunks { done: Some(1) }
                }
            }
        );

        saves.one_chunk_done(2);
        assert_eq!(saves.inner, Chunks { done: Some(2) });
        assert_eq!(
            Saves::<Chunks>::detect_crash(save_dir, 3),
            SavesState::Chunk {
                saves: Saves {
                    path: PathBuf::from(save_dir),
                    inner: Chunks { done: Some(2) }
                }
            }
        );

        let mut saves = saves.chunks_all_done();
        assert_eq!(saves.inner.merges_done, HashSet::new());
        assert_eq!(
            Saves::<Chunks>::detect_crash(save_dir, 3),
            SavesState::Merge {
                saves: Saves {
                    path: PathBuf::from(save_dir),
                    inner: Merge {
                        merges_done: HashSet::new()
                    }
                }
            }
        );

        // simulates a few merge
        let first_merge_done = [8, 6, 9];

        for i in first_merge_done {
            saves.one_merge_done(i);
        }
        assert_eq!(
            Saves::<Chunks>::detect_crash(save_dir, 3),
            SavesState::Merge {
                saves: Saves {
                    path: PathBuf::from(save_dir),
                    inner: Merge {
                        merges_done: HashSet::from_iter(first_merge_done)
                    }
                }
            }
        );

        for i in (0..NB_FILE_IN_AN_INDEX).filter(|x| !first_merge_done.contains(x)) {
            saves.one_merge_done(i);
        }
        saves.merge_all_done();
    }

    #[rstest]
    fn test_store_and_load_atomics(random_directory: AutoRemoveDirectory) {
        let save_dir = random_directory.filename().to_str().unwrap();
        std::fs::create_dir_all(save_dir).expect("Failed to create test directory");
        let save_dir = Path::new(save_dir);

        let atomic_total_kmers = AtomicU64::new(85);
        let atomic_dense_kmers_count = AtomicU64::new(89485);
        let atomic_sparse_kmers_count = AtomicU64::new(7523);

        let mut saves = Saves::<Chunks>::new(save_dir.to_path_buf());

        for i in 0..8 {
            saves.one_chunk_done(i);
        }

        assert_eq!(saves.inner, Chunks { done: Some(7) });

        saves.store_current_atomics_to_disk(
            &atomic_total_kmers,
            &atomic_dense_kmers_count,
            &atomic_sparse_kmers_count,
        );
        saves.one_chunk_done(8);
        let (a, b, c) = saves.load_atomics_from_disk();

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

    #[rstest]
    fn test_store_and_load_debug_atomics(random_directory: AutoRemoveDirectory) {
        let save_dir = random_directory.filename().to_str().unwrap();
        std::fs::create_dir_all(save_dir).expect("Failed to create test directory");
        let save_dir = Path::new(save_dir);

        let mut saves = Saves::<Chunks>::new(save_dir.to_path_buf());

        for i in 0..8 {
            saves.one_chunk_done(i);
        }

        assert_eq!(saves.inner, Chunks { done: Some(7) });

        let atomic_sparse_one_seen = std::sync::atomic::AtomicU64::new(85);
        let atomic_sparse_fp_seen = std::sync::atomic::AtomicU64::new(89485);

        saves.store_current_chunk_debug_atomics_to_disk(
            &atomic_sparse_one_seen,
            &atomic_sparse_fp_seen,
        );
        saves.one_chunk_done(8);
        let (a, b) = saves.load_debug_atomics_from_disk();

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
