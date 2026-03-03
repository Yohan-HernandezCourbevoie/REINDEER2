/// Returns the minimum number of partitions required to index `nb_files` files with `nb_abundance_level` abundance levels and filters of size `bf_size`.
#[must_use]
pub fn get_number_of_partitions(nb_files: usize, nb_abundance_level: usize, bf_size: u64) -> usize {
    let max_size_roaring = 2u64.pow(32);

    usize::try_from(
        (nb_abundance_level as u64 * nb_files as u64 * bf_size).div_ceil(max_size_roaring),
    )
    .expect("error in conversion")
}

/// If `requested_number_of_files` is enough to index files `files_paths`, returns `requested_number_of_files`. Returns `file_paths.len()` otherwise.
#[must_use]
pub fn get_min_number_of_files<T>(
    file_paths: &[T],
    requested_number_of_files: Option<usize>,
) -> usize {
    let number_of_files = file_paths.len();
    match requested_number_of_files {
        Some(requested_number_of_files) => {
            if requested_number_of_files < number_of_files {
                let msg = format!(
                "The requested number of file was {}, but there are {} files to index. Continuing with space for exactly {} files. \
                If you want to be able to add more datasets, increase the number of requested files.",
                requested_number_of_files, number_of_files, number_of_files
            );
                log::warn!("{msg}");
                println!("{msg}");
                number_of_files
            } else {
                requested_number_of_files
            }
        }
        None => number_of_files,
    }
}
