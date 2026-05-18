use itertools::Itertools;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::{
    cmp::min,
    fs::{self, File},
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::Path,
};

use crate::reindeer2::{
    NB_FILE_IN_AN_INDEX, Parameters, Reindeer2, create_and_reserve_tar_get_file,
    merge::merge_partition_slices_interleaved,
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

    assert!(
        new_color_number <= index_ref.parameters.capacity,
        "not enough capacity to store all the indexes"
    );

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

    #[allow(clippy::map_flatten, reason = "readability")]
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
