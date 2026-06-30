use itertools::Itertools;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::{
    cmp::min,
    fs::File,
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
    sync::{Arc, Mutex},
    time::Instant,
};

#[cfg(feature = "self-destruct")]
use crate::FailIndexation;

use crate::reindeer2::{
    NB_FILE_IN_AN_INDEX, create_and_reserve_tar_get_file,
    merge::merge_partition_slices_interleaved,
    saves::{Merge, Saves},
    storage::filters::load_bloom_filter_from_big_file,
};
/// Removes partitions (they are unecessary after a merge).
pub fn remove_merged_partitions(
    chunk_files_dir: &str,
    nb_chunk: usize, // number of colors in each chunk
) -> io::Result<()> {
    (0..NB_FILE_IN_AN_INDEX)
        .into_par_iter()
        .try_for_each(|file_id| {
            for chunk_id in 0..nb_chunk {
                let path = format!(
                    "{chunk_files_dir}/partition_bloom_filters_group{file_id}_chunk{chunk_id}.bin"
                );
                std::fs::remove_file(path)?;
            }
            Ok::<(), io::Error>(())
        })?;
    Ok(())
}

/// Merges all partitions of an index being built.
pub fn merge_all_partitions_of_chunks(
    chunk_files_dir: &Path,
    output_dir: &Path,
    partitioned_bf_size: usize,
    abundance_number: usize,
    color_counts_per_chunk: &[usize], // number of colors in each chunk
    num_partitions: usize,
    mut saves: Saves<Merge>,
    #[cfg(feature = "self-destruct")] fail: &Option<FailIndexation>,
) -> io::Result<()> {
    let start_time = Instant::now();

    let nb_chunk = color_counts_per_chunk.len();
    let nb_partition_in_a_file = num_partitions.div_ceil(NB_FILE_IN_AN_INDEX);
    let still_to_be_merge = saves.get_still_to_be_merged();
    let saves_arc = Arc::new(Mutex::new(&mut saves));

    still_to_be_merge.into_par_iter().try_for_each(|file_id| {
        let range_start = nb_partition_in_a_file * file_id;
        let range_end = min(nb_partition_in_a_file * (file_id + 1), num_partitions);
        let mut range = range_start..range_end;

        let path =
            Path::new(output_dir).join(format!("partition_bloom_filters_group{file_id}.bin"));
        let len: usize = range
            .try_len()
            .expect("range object should have a known length");
        let len =
            u64::try_from(len).expect("should have less than u64::MAX partitions in a single file");
        let mut file = create_and_reserve_tar_get_file(&path, len);

        let inpaths = (0..nb_chunk)
            .map(|chunk_id| {
                chunk_files_dir.join(format!(
                    "partition_bloom_filters_group{file_id}_chunk{chunk_id}.bin"
                ))
            })
            .collect_vec();

        let merge_result = range.try_for_each(|partition_idx| {
            let index = (partition_idx % nb_partition_in_a_file) as u64;
            merge_into_single_tar_get_file(
                &inpaths,
                partitioned_bf_size,
                abundance_number,
                color_counts_per_chunk,
                index,
                &mut file,
            )?;
            Ok::<(), io::Error>(())
        });

        let mut saves = saves_arc.lock().expect(
            "fatal error: a thread holding the mutex panicked, so this thread will panic as well",
        );
        #[cfg(feature = "self-destruct")]
        if let Some(FailIndexation::Merge(fail_merge)) = fail
            && *fail_merge == file_id
        {
            panic!("indexation failed on merge {fail_merge} as planned")
        }
        saves.one_merge_done(file_id);

        // remove the intermedite partitions after the merge
        for path in inpaths {
            std::fs::remove_file(path)?;
        }

        merge_result
    })?;

    let elapsed_time = start_time.elapsed();
    log::info!(
        "All partitions merged and written to disk in {:.2?}",
        elapsed_time
    );
    saves.merge_all_done();

    Ok(())
}

/// Merges all filter from tar_get files chunk_files (on tar_get file per chunk) into a unique tar_get file.
fn merge_into_single_tar_get_file(
    chunk_files: &[PathBuf],
    partitioned_bf_size: usize,
    abundance_number: usize,
    color_counts: &[usize], // number of colors for each chunk
    index: u64,
    output_file: &mut File,
) -> io::Result<RoaringBitmap> {
    if chunk_files.is_empty() || chunk_files.len() != color_counts.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Mismatch between chunk files and color counts, or no chunks provided.",
        ));
    }

    // OPTIMIZE we may be able to load only one filter at a time
    let filters: Vec<(RoaringBitmap, usize)> = chunk_files
        .iter()
        .zip(color_counts)
        .map(
            |(chunk_file, nb_color)| match load_bloom_filter_from_big_file(chunk_file, index) {
                Ok(filter) => (filter, *nb_color),
                Err(e) => {
                    log::error!(
                        "Failed to load Bloom filter {}: {}",
                        chunk_file.display(),
                        e
                    );
                    panic!(
                        "Failed to load Bloom filter {}: {}",
                        chunk_file.display(),
                        e
                    );
                }
            },
        )
        .collect();

    // write all slices in the right order in a  larger filter
    let final_bf =
        merge_partition_slices_interleaved(&filters, partitioned_bf_size, abundance_number);

    // write as a targ_get file
    let serializer =
        |writer: &mut BufWriter<&mut File>, bitmap: &RoaringBitmap| bitmap.serialize_into(writer);
    // let mut output_file = output_file.lock().unwrap();
    tar_get::append_element(output_file, &final_bf, serializer).map_err(|err| {
        io::Error::other(format!("Error while eppending to a tar_get file: {err}"))
    })?;
    output_file.flush()?;
    Ok(final_bf)
}

#[cfg(test)]
mod tests {

    use std::fs::OpenOptions;

    use super::*;

    use crate::reindeer2::test_utils::AutoRemoveDirectory;
    use rstest::{fixture, rstest};

    #[fixture]
    pub fn random_directory() -> AutoRemoveDirectory {
        AutoRemoveDirectory::create_random()
    }

    #[rstest]
    fn test_merge_partition_bloom_filters(random_directory: AutoRemoveDirectory) {
        use roaring::RoaringBitmap;
        use std::fs::create_dir_all;

        let test_dir = random_directory.filename();
        create_dir_all(test_dir).expect("Failed to create test directory");

        let chunk1_path = test_dir.join("chunk1_p0.bin");
        let chunk2_path = test_dir.join("chunk2_p0.bin");
        let chunk3_path = test_dir.join("chunk3_p0.bin");

        let partitioned_bf_size = 2;
        let abundance_number = 3;
        // test Bloom filters
        let mut bf1 = RoaringBitmap::new();
        bf1.insert(1);
        bf1.insert(2); //011000 000000
        let mut bf2 = RoaringBitmap::new();
        bf2.insert(3);
        bf2.insert(4); // 000110 000000
        let mut bf3 = RoaringBitmap::new();
        bf3.insert(5); // 000 001

        //expected
        // 011000000110000 000000000000001 [1,2,9,10,29]

        let chunk_files = vec![
            chunk1_path.clone(),
            chunk2_path.clone(),
            chunk3_path.clone(),
        ];
        let color_counts = vec![2, 2, 1]; //nb of colors for each chunk
        // only on file per tar_get file
        let len = 1;

        for (chunk_path, (bf, _colors)) in
            chunk_files.iter().zip(vec![(bf1, 2), (bf2, 2), (bf3, 1)])
        {
            let mut file = create_and_reserve_tar_get_file(chunk_path, len);
            tar_get::append_element(&mut file, &bf, |writer, bitmap| {
                bitmap.serialize_into(writer)
            })
            .expect("should have been able to add an element to a tar_get file");
        }

        let path = Path::new(test_dir).join("partition_bloom_filters.bin");
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&path)
            .unwrap();
        let mut writer = BufWriter::new(file);
        tar_get::reserve_capacity(&mut writer, len).unwrap(); // TODO expect

        let mut file = writer.into_inner().unwrap();

        // Call the modified function
        let merged_bf = merge_into_single_tar_get_file(
            &chunk_files,
            partitioned_bf_size,
            abundance_number,
            &color_counts,
            0, // since there is only one file par tar_get file, we want the one at index 0
            &mut file,
        )
        .expect("Failed to merge partition Bloom filters");

        // check the merged Bloom filter
        let expected_elements: Vec<u32> = vec![1, 2, 9, 10, 29];
        for elem in expected_elements {
            assert!(
                merged_bf.contains(elem),
                "Expected element {} not found in the merged Bloom filter",
                elem
            );
        }
    }
}
