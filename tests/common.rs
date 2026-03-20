//! Common utilities for integration tests

use std::io::BufRead;

use uuid::Uuid;

/// Checks is all the elements of a slice are the same.
pub fn is_all_same<T: PartialEq>(arr: &[T]) -> bool {
    if arr.is_empty() {
        return true;
    }
    let first = &arr[0];
    arr.iter().all(|item| item == first)
}

fn count_lines<P: AsRef<std::path::Path>>(path: P) -> std::io::Result<usize> {
    let mut count = 0;
    for line in std::io::BufReader::new(std::fs::File::open(path)?).lines() {
        line?;
        count += 1;
    }
    Ok(count)
}

/// Checks if the given number of partitions would lead to a crash at indexation time.
/// Returns the number of partitions, or, in case of a crash, returns the minimum number of partition that is safe to use.
#[must_use]
pub fn check_number_of_partitions<P: AsRef<std::path::Path>>(
    fof: P,
    nb_partitions: usize,
    nb_abundance_level: usize,
    bf_size: u64,
) -> usize {
    let nb_files = count_lines(&fof).unwrap();
    let nb_min_partitions = nb_abundance_level
        * nb_files
        * usize::try_from(bf_size / 2u64.pow(32)).expect("error in conversion");

    if nb_partitions < nb_min_partitions {
        nb_min_partitions
    } else {
        nb_partitions
    }
}

/// Wraps a directory and remove it when it is `drop`.
pub struct AutoRemoveDirectory {
    path: String,
}

impl AutoRemoveDirectory {
    /// Creates a new random direcory that will be removed on `drop`.
    pub fn create_random() -> Self {
        let my_uuid = Uuid::new_v4();
        let filename = format!("{my_uuid}");
        Self { path: filename }
    }

    // /// Wraps a direcory that will be removed on `drop`.
    // /// # Warning
    // /// This will remove the directory.
    // pub fn create_from_path(path: String) -> Self {
    //     Self { path }
    // }

    pub fn filename(&self) -> &str {
        &self.path
    }
}

impl Drop for AutoRemoveDirectory {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.path).unwrap()
    }
}

/// Wraps a file and remove it when it is `drop`.
pub struct AutoRemoveFile {
    path: String,
}

impl AutoRemoveFile {
    // /// Creates a new random file that will be removed on `drop`.
    // pub fn create_random() -> Self {
    //     let my_uuid = Uuid::new_v4();
    //     let filename = format!("{my_uuid}");
    //     Self { path: filename }
    // }

    /// Wraps a file that will be removed on `drop`.
    /// # Warning
    /// This will remove the file.
    pub fn create_from_path(path: String) -> Self {
        Self { path }
    }

    pub fn filename(&self) -> &str {
        &self.path
    }
}

impl Drop for AutoRemoveFile {
    fn drop(&mut self) {
        std::fs::remove_file(&self.path).unwrap()
    }
}
