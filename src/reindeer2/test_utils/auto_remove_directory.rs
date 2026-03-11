use std::path::{Path, PathBuf};

/// Wraps a directory and remove it when it is `drop`.
#[derive(Debug)]
pub struct AutoRemoveDirectory {
    path: PathBuf,
}

impl AutoRemoveDirectory {
    /// Creates a new random direcory name that will be removed on `drop`.
    /// The directory is *not* created, but is expected to exist before `drop`.
    pub fn create_random() -> Self {
        let my_uuid = uuid::Uuid::new_v4();
        let filename = format!("{my_uuid}");
        Self {
            path: PathBuf::from(filename),
        }
    }

    // /// Wraps a direcory that will be removed on `drop`.
    // /// # Warning
    // /// This will remove the directory.
    // pub fn create_from_path(path: String) -> Self {
    //     Self { path }
    // }

    /// Gets the filename of the directory. Warning: any object file placed inside will be removed on `drop`.
    pub fn filename(&self) -> &Path {
        &self.path
    }
}

impl Drop for AutoRemoveDirectory {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.path).unwrap()
    }
}
