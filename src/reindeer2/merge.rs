use itertools::Itertools;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::{
    cmp::min,
    fs::{self, File},
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::Path,
    time::Instant,
};

use crate::reindeer2::{
    create_and_reserve_tar_get_file, filter::load_raw_bloom_filter, Parameters, Reindeer2,
    NB_FILE_IN_AN_INDEX,
};

fn read_lines_of_file_of_indexes(indexes_fof: &str) -> io::Result<Vec<String>> {
    let file = File::open(indexes_fof)?;
    Ok(BufReader::new(file)
                .lines()
                .map(|line| {
                    line.unwrap_or_else(|_| {
                        let msg = format!("Fatal error: when merging indexes, an error was encountered when reading the file {indexes_fof}.");
                        println!("{msg}");  // warn user
                        log::error!("{msg}");  // log error
                        // TODO discuss if we panic or exit
                        panic!("should have been able to read file {}", indexes_fof);  // panic
                    })
                })
                .map(|line| line.trim().to_string())
                .filter(|line| !line.is_empty())
                .collect())
}

/// Merges `tar_get` files `tar_get_file_id` from all indexes in `indexes` and produces a single tar_get file in the `output_dir` directory.
fn merge_tar_get_files(
    indexes: &[Reindeer2],
    output_dir: &str,
    tar_get_file_id: usize,
    nb_files_in_tar_get: u64,
    partitioned_bf_size: usize,
    abundance_number: usize,
) {
    let tar_get_file_name = format!("partition_bloom_filters_group{tar_get_file_id}.bin");
    let mut tar_get_readers_and_colors = indexes
        .iter()
        .map(|index| {
            let path = Path::new(&index.index_dir).join(&tar_get_file_name);
            let tar_get_file = File::open(&path).unwrap_or_else(|_| {
                panic!("should have been able to open file {}", path.display())
            });
            (BufReader::new(tar_get_file), index.parameters.nb_color)
        })
        .collect_vec();

    // FIXME check the numer of file inside all tar_get files

    let tar_get_path = Path::new(&output_dir).join(&tar_get_file_name);
    let mut tar_get_writer = create_and_reserve_tar_get_file(&tar_get_path, nb_files_in_tar_get);

    let deserializer = |reader: Vec<u8>| {
        let slice: &[u8] = &reader;
        RoaringBitmap::deserialize_from(slice)
    };

    let serializer =
        |writer: &mut BufWriter<&mut File>, bitmap: &RoaringBitmap| bitmap.serialize_into(writer);

    for index in 0..nb_files_in_tar_get {
        // extract all files i from all tar_get files
        // TODO error handling
        let filters = tar_get_readers_and_colors
            .iter_mut()
            .enumerate()
            .map(|(i, (reader, nb_color))| {
                let filter =
                    tar_get::deserialize(reader, index, deserializer).unwrap_or_else(|err| {
                        panic!(
                            "should have been able to deserialize file {} from tar_get file {} from index {} ({})",
                            index, tar_get_file_id, indexes[i].index_dir, err
                        )
                    });
                (filter, *nb_color)
            })
            .collect_vec();

        // write all slices in the right order in a  larger filter
        let final_bf = merge_partition_slices_interleaved(
            &filters,
            partitioned_bf_size,
            abundance_number,
            // color_counts,
        );

        // write filter to disk
        tar_get::append_element(&mut tar_get_writer, &final_bf, serializer).unwrap_or_else(|_| {
            panic!(
                "shoudl have been able to flush file {}",
                tar_get_path.display()
            )
        });
        tar_get_writer.flush().unwrap_or_else(|_| {
            panic!(
                "shoudl have been able to flush file {}",
                tar_get_path.display()
            )
        });
    }
}

/// Merges an arbitrary number of indexes listed in the file `indexes_fof`. The new index is placed in `output_dir`.
///
/// Each line of the fof is expected to the path to one index directory.
pub fn merge_multiple_indexes(indexes_fof: &str, output_dir: &str) -> io::Result<()> {
    // read the list of index directories.
    let index_dirs = read_lines_of_file_of_indexes(indexes_fof)?;

    if index_dirs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No index directories provided in the file-of-indexes",
        ));
    }
    if index_dirs.len() == 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Only one index directory provided, nothing to do",
        ));
    }

    let index_ref = Reindeer2::load_from_disk(&index_dirs[0]).unwrap_or_else(|_| {
        panic!(
            "should have been able to load index infos from disk {}",
            &index_dirs[0]
        )
    });

    //  vector to store (index_dir, color_count) for every index.
    let indexes = index_dirs
        .iter()
        .map(|index_dir| {
            Reindeer2::load_from_disk(index_dir).unwrap_or_else(|_| {
                panic!(
                    "should have been able to load index infos from disk {}",
                    index_dir
                )
            })
        })
        .collect_vec();

    // check that the sum of all nb_files is smaller than the capacity

    // let mut indexes_metadata = vec![(index_dirs[0].clone(), index_ref.parameters.nb_color)];
    let mut new_color_number = index_ref.parameters.nb_color;
    assert!(
        new_color_number <= index_ref.parameters.capacity,
        "not enough capacity to store all the indexes"
    );

    // for all other indexes, check that the parameters match and add its color count
    for index_to_merge in indexes.iter().skip(1) {
        if !index_ref.parameters.can_merge(&index_to_merge.parameters) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Index {} does not match parameters of the first index ({})",
                    index_to_merge.index_dir, index_ref.index_dir
                ),
            ));
        }
        new_color_number += index_to_merge.parameters.nb_color;
    }

    fs::create_dir_all(output_dir)?;

    let partitioned_bf_size =
        (index_ref.parameters.bf_size as usize) / index_ref.parameters.partition_number;

    let nb_partition = index_ref.parameters.partition_number;
    let abundance_number: usize = index_ref.parameters.abundance_number.into();

    let nb_partition_in_a_file = nb_partition.div_ceil(NB_FILE_IN_AN_INDEX);

    // too unlikely for me to bother right now
    let nb_partition_in_a_file = u64::try_from(nb_partition_in_a_file).expect("Number of partition is to big. Should have less than u64::MAX partitions in a single file. Please contact the authors to allow for more partitions.");

    let nb_partition = u64::try_from(nb_partition).expect("Number of partition is to big. Should have less than u64::MAX partitions. Please contact the authors to allow for more partitions.");

    (0..NB_FILE_IN_AN_INDEX)
        .into_par_iter()
        .for_each(|tar_get_file_id| {
            let tar_get_file_id_u64 = u64::try_from(tar_get_file_id)
                .expect("Too many tar_get_files in a single index. Please contact the authors.");

            let range_start = nb_partition_in_a_file * tar_get_file_id_u64;
            let range_end = min(
                nb_partition_in_a_file * (tar_get_file_id_u64 + 1),
                nb_partition,
            );

            if range_end > range_start {
                let nb_files_in_tar_get = range_end - range_start;

                merge_tar_get_files(
                    &indexes,
                    output_dir,
                    tar_get_file_id,
                    nb_files_in_tar_get,
                    partitioned_bf_size,
                    abundance_number,
                );
            }
        });

    let parameters = Parameters {
        nb_color: new_color_number,
        ..index_ref.parameters
    };

    let indexed_file_names = indexes
        .into_iter()
        .map(|index| index.indexed_file_names)
        .flatten()
        .collect_vec();

    let index_merged = Reindeer2 {
        parameters,
        indexed_file_names,
        index_dir: String::from(output_dir),
    };
    index_merged.save_to_disk()?;

    log::info!(
        "Successfully merged {} indexes into directory: {}",
        index_dirs.len(),
        output_dir
    );
    Ok(())
}

// // --- MERGE FUNCTIONS ---

/// Merges all partitions of an index being built.
///
/// The initial partitions
pub fn merge_all_partitions_of_chunks(
    chunk_files_dir: &str,
    output_dir: &str,
    partitioned_bf_size: usize,
    abundance_number: usize,
    color_counts_per_chunk: Vec<usize>, // number of colors in each chunk
    num_partitions: usize,
) -> io::Result<()> {
    let start_time = Instant::now();

    let nb_chunk = color_counts_per_chunk.len();
    let nb_partition_in_a_file = num_partitions.div_ceil(NB_FILE_IN_AN_INDEX);

    (0..NB_FILE_IN_AN_INDEX)
        .into_iter()
        .try_for_each(|file_id| {
            let range_start = nb_partition_in_a_file * file_id;
            let range_end = min(nb_partition_in_a_file * (file_id + 1), num_partitions);
            let mut range = range_start..range_end;

            let path =
                Path::new(output_dir).join(format!("partition_bloom_filters_group{file_id}.bin"));
            let len: usize = range
                .try_len()
                .expect("range object should have a known length");
            let len = u64::try_from(len)
                .expect("should have less than u64::MAX partitions in a single file");
            let mut file = create_and_reserve_tar_get_file(&path, len);

            range.try_for_each(|partition_idx| {
                // collect chunk files for the current partition
                let chunk_files_for_partition: Vec<String> = (0..nb_chunk)
                    .map(|chunk_idx| {
                        format!(
                            "{}/partition_bloom_filters_c{}_p{}.bin",
                            chunk_files_dir, chunk_idx, partition_idx
                        )
                    })
                    .collect();

                // merge all chunks for the current partition + serialize
                merge_partition_bloom_filters(
                    &chunk_files_for_partition,
                    partitioned_bf_size,
                    abundance_number,
                    &color_counts_per_chunk,
                    &mut file,
                )?;

                Ok::<(), io::Error>(())
            })
        })?;

    let elapsed_time = start_time.elapsed();
    log::info!(
        "All partitions merged and written to disk in {:.2?}",
        elapsed_time
    );

    Ok(())
}

/// Merges all filter from files chunk_files into a unique tar_get file.
pub fn merge_partition_bloom_filters(
    chunk_files: &[String],
    partitioned_bf_size: usize,
    abundance_number: usize,
    color_counts: &[usize], // number of colors for each chunk
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
        .map(|chunk_file| match load_raw_bloom_filter(chunk_file) {
            Ok(filter) => filter,
            Err(e) => {
                log::error!("Failed to load Bloom filter {}: {}", chunk_file, e);
                panic!("Failed to load Bloom filter {}: {}", chunk_file, e);
            }
        })
        .collect();

    // write all slices in the right order in a  larger filter
    let final_bf = merge_partition_slices_interleaved(
        &filters,
        partitioned_bf_size,
        abundance_number,
        // color_counts,
    );

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

/* example for build_new_bitset_with_gaps_from_merged + interleave_slices_with_zero_runs
 merged = 110001 000000 (2 "colums" in the partition, 2 datasets, 3 abundances)
 bf to add = 111 000 (2 "colums" in the partition, 1 dataset, 3 abundances)

 build_new_bitset_with_gaps_from_merge will prepare merge as follows:

 110001 000 000000 000 <- new runs of 0 of size 3 to add abundances for the new dataset brought by bf_to_add

 conversely, interleave_slices_with_zero_runs will prepare bf_to_add for the union:

 000000 111 000000 000 <- new runs of 0 of size 2*3 (size of a block in the merged vector)

 then we perform the union

 110001111 000000000 -> 2 "colums" in the partition, THREE datasets, 3 abundances

*/

/// Merges all filter from files `chunk_files` into a `RoaringBitmap`.
fn merge_partition_slices_interleaved(
    filters: &[(RoaringBitmap, usize)],
    partitioned_bf_size: usize,
    abundance_number: usize,
    // color_counts: &[usize], // number of colors in each chunk TODO enlever
) -> RoaringBitmap {
    // preload all Bfs

    // vector to collect all positions for the final bf
    let mut final_positions = Vec::new();

    let sum_color: usize = filters.iter().map(|(_filter, nb_color)| nb_color).sum();

    for slice_idx in 0..partitioned_bf_size {
        let mut current_offset = slice_idx * abundance_number * sum_color;

        for (filter, nb_color) in filters {
            let slice_start = slice_idx * abundance_number * nb_color;
            let slice_end = slice_start + abundance_number * nb_color;

            let slice_start_u32 = slice_start as u32;
            let current_offset_u32 = current_offset as u32;
            let slice_end_u32 = slice_end as u32;
            #[cfg(any(debug_assertions, test))]
            {
                assert_eq!(slice_start, slice_start_u32 as usize);
                assert_eq!(current_offset, current_offset_u32 as usize);
                assert_eq!(slice_end, slice_end_u32 as usize);
            }
            // collect positions in the slice and adjust by offset
            final_positions.extend(
                filter
                    .range(slice_start_u32..slice_end_u32)
                    .map(|pos| pos - slice_start_u32 + current_offset_u32),
            );

            // update the offset for the next chunk in this slice
            current_offset += abundance_number * nb_color;
        }
    }
    RoaringBitmap::from_sorted_iter(final_positions)
        .expect("Attempt to merge with unsorted positions")
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
        use std::fs::{create_dir_all, File};
        use std::io::Write;

        let test_dir = random_directory.filename().to_str().unwrap();
        create_dir_all(test_dir).expect("Failed to create test directory");

        let chunk1_path = format!("{}/chunk1_p0.bin", test_dir);
        let chunk2_path = format!("{}/chunk2_p0.bin", test_dir);
        let chunk3_path = format!("{}/chunk3_p0.bin", test_dir);

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

        for (chunk_path, (bf, colors)) in chunk_files.iter().zip(vec![(bf1, 2), (bf2, 2), (bf3, 1)])
        {
            let mut file = File::create(chunk_path).expect("Failed to create test chunk file");
            file.write_all(&(colors as u64).to_le_bytes())
                .expect("Failed to write color count");
            bf.serialize_into(&mut file)
                .expect("Failed to serialize Bloom filter");
        }

        let path = Path::new(test_dir).join("partition_bloom_filters.bin");
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&path)
            .unwrap();
        let mut writer = BufWriter::new(file);
        let len = 1;
        tar_get::reserve_capacity(&mut writer, len).unwrap(); // TODO expect

        let mut file = writer.into_inner().unwrap();

        // Call the modified function
        let merged_bf = merge_partition_bloom_filters(
            &chunk_files,
            partitioned_bf_size,
            abundance_number,
            &color_counts,
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
