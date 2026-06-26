//! Saves states for REINDEER2

use std::{
    collections::HashSet,
    path::{Path, PathBuf},
};

use crate::reindeer2::NB_FILE_IN_AN_INDEX;

/// Saves in the middle of doing chunks
#[derive(Debug, PartialEq, Eq)]
pub struct Chunks {
    done: Option<usize>,
}

/// Saves in the middle of doing chunks
pub struct Merge {
    inner: HashSet<usize>,
}

/// Saves of the index
pub struct Saves<T> {
    path: PathBuf,
    inner: T,
}

#[derive(Debug, PartialEq)]
pub enum CrashState {
    Chunk { last_done: Option<usize> },
    Merge { still_to_be_done: HashSet<usize> },
}

#[derive(Copy, Clone)]
pub struct CrashedChunk {
    pub last_done: Option<usize>,
}

#[derive(Clone)]
pub struct CrashedMerge {
    pub still_to_be_done: HashSet<usize>,
}

impl<T> Saves<T> {
    pub fn get_filename_for_chunk(chunk: usize) -> String {
        format!("chunk_{chunk}_done")
    }

    pub fn get_filename_for_all_chunk_done() -> &'static str {
        "chunks_all_done"
    }

    pub fn get_filename_for_merge(merge: usize) -> String {
        format!("merge_{merge}_done")
    }

    pub fn detect_crash(path: &Path, nb_chunk: usize) -> CrashState {
        let crashed_after_chunk =
            std::fs::exists(path.join(Self::get_filename_for_all_chunk_done()))
                .expect("should have been able to read the save folder");
        if crashed_after_chunk {
            let mut merge_to_be_done = HashSet::new();
            for i in 0..NB_FILE_IN_AN_INDEX {
                let chunk_filename = Self::get_filename_for_merge(i);
                let chunk_path = path.join(&chunk_filename);
                let chunk_done = std::fs::exists(chunk_path)
                    .expect("should have been able to check if the chunk was done");
                if !chunk_done {
                    merge_to_be_done.insert(i);
                }
            }

            assert!(!merge_to_be_done.is_empty());

            CrashState::Merge {
                still_to_be_done: merge_to_be_done,
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

            CrashState::Chunk { last_done }
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

    pub const fn from_crashed_chunk(path: PathBuf, crash_state: Option<CrashedChunk>) -> Self {
        match crash_state {
            None => Self::new(path),
            Some(CrashedChunk { last_done }) => Self {
                path,
                inner: Chunks { done: last_done },
            },
        }
    }

    pub const fn get_nb_chunk_that_can_be_skipped(&self) -> usize {
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
                inner: HashSet::new(),
            },
        }
    }
}

impl Saves<Merge> {
    // panics if the merge was already done
    pub fn one_merge_done(&mut self, merge_done: usize) {
        let inserted = self.inner.inner.insert(merge_done);
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
            CrashState::Chunk { last_done: None }
        );

        saves.one_chunk_done(0);
        assert_eq!(saves.inner, Chunks { done: Some(0) });
        assert_eq!(
            Saves::<Chunks>::detect_crash(save_dir, 3),
            CrashState::Chunk { last_done: Some(0) }
        );

        saves.one_chunk_done(1);
        assert_eq!(saves.inner, Chunks { done: Some(1) });
        assert_eq!(
            Saves::<Chunks>::detect_crash(save_dir, 3),
            CrashState::Chunk { last_done: Some(1) }
        );

        saves.one_chunk_done(2);
        assert_eq!(saves.inner, Chunks { done: Some(2) });
        assert_eq!(
            Saves::<Chunks>::detect_crash(save_dir, 3),
            CrashState::Chunk { last_done: Some(2) }
        );

        let mut saves = saves.chunks_all_done();
        assert_eq!(saves.inner.inner, HashSet::new());
        assert_eq!(
            Saves::<Chunks>::detect_crash(save_dir, 3),
            CrashState::Merge {
                still_to_be_done: HashSet::from_iter(0..NB_FILE_IN_AN_INDEX)
            }
        );

        // simulates a few merge
        let first_merge_done = [8, 6, 9];
        let still_to_be_done = (0..NB_FILE_IN_AN_INDEX)
            .filter(|x| !first_merge_done.contains(x))
            .collect();
        for i in first_merge_done {
            saves.one_merge_done(i);
        }
        assert_eq!(
            Saves::<Chunks>::detect_crash(save_dir, 3),
            CrashState::Merge { still_to_be_done }
        );

        for i in (0..NB_FILE_IN_AN_INDEX).filter(|x| !first_merge_done.contains(x)) {
            saves.one_merge_done(i);
        }
        saves.merge_all_done();
    }
}
