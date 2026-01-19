use colored::Colorize;
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};
fn count_lines<P: AsRef<Path>>(path: P) -> std::io::Result<usize> {
    let mut count = 0;
    for line in BufReader::new(File::open(path)?).lines() {
        line?;
        count += 1;
    }
    Ok(count)
}

/// Checks if the given number of partitions would lead to a crash at indexation time.
/// Returns the number of partitions, or, in case of a crash, returns the minimum number of partition that is safe to use.
#[must_use]
pub fn check_number_of_partitions<P: AsRef<Path>>(
    fof: P,
    nb_partitions: usize,
    nb_abundance_level: usize,
    bf_size: u64,
) -> usize {
    let nb_files = count_lines(&fof).unwrap_or_else(|_| {
        panic!(
            "should have been able to read {}",
            fof.as_ref().to_str().unwrap()
        )
    });
    let nb_min_partitions = nb_abundance_level
        * nb_files
        * usize::try_from(bf_size / 2u64.pow(32)).expect("error in conversion");

    if nb_partitions < nb_min_partitions {
        log::warn!(
            "{}",
            format!(
                "Warning: the requested number of partition was {}. \
                Keeping this value would lead to a crash. \
                The minimum number of partitions is {}. \
                Indexation will continue with this minimum number of partitions to prevent crash. \
                Using the minimum number of partitions prevents adding datasets in the index in the future. \
                If you want to be able to add more datasets, increase the number of partitions.",
                nb_partitions, nb_min_partitions
            )
            .yellow()
        );
        nb_min_partitions
    } else {
        nb_partitions
    }
}
