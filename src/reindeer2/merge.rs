use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::{
    fs::{self, File, OpenOptions},
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::Path,
    sync::{Arc, Mutex},
    time::Instant,
};

use crate::reindeer2::{filter::load_raw_bloom_filter, Parameters, Reindeer2};

// merge an arbitrary number of indexes  listed in a fof
// each line of the fof is expected to he path to one index dir
// in the end, alsos write a  new metadata JSON file in the output dir
pub fn merge_multiple_indexes(indexes_fof: &str, output_dir: &str) -> io::Result<()> {
    // read the list of index directories.
    let index_dirs: Vec<String> = {
        let file = File::open(indexes_fof)?;
        BufReader::new(file)
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
            .collect()
    };

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

    // read metadata from the first index as base parameter
    // let (k, m, bf_size, partition_number, first_color, abundance_number, abundance_max, dense_option) =
    // read_partition_from_csv(&index_dirs[0], "index_info.csv")?;
    let index_ref = Reindeer2::load_from_disk(&index_dirs[0]).expect(&format!(
        "should have been able to load index infos from disk {}",
        &index_dirs[0]
    ));

    //  vector to store (index_dir, color_count) for every index.
    let mut indexes_metadata = vec![(index_dirs[0].clone(), index_ref.parameters.nb_color)];
    let mut new_color_number = index_ref.parameters.nb_color;
    let mut indexed_file_names = index_ref.indexed_file_names;

    // for all other indexes, check that the parameters match and add its color count
    for index_dir in index_dirs.iter().skip(1) {
        let index_to_merge = Reindeer2::load_from_disk(index_dir).expect(&format!(
            "should have been able to load index infos from disk {}",
            index_dir
        ));
        if !index_ref.parameters.can_merge(&index_to_merge.parameters) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Index {} does not match parameters of the first index ({})",
                    index_dir, index_ref.index_dir
                ),
            ));
        }
        new_color_number += index_to_merge.parameters.nb_color;
        indexes_metadata.push((index_dir.clone(), index_to_merge.parameters.nb_color));
        indexed_file_names.extend(index_to_merge.indexed_file_names);
    }

    fs::create_dir_all(output_dir)?;

    let partitioned_bf_size =
        (index_ref.parameters.bf_size as usize) / index_ref.parameters.partition_number;

    let path = Path::new(output_dir).join("partition_bloom_filters.bin");
    let file = OpenOptions::new()
        .read(true)
        .write(true)
        .create_new(true)
        .open(&path)
        .unwrap();
    let mut writer = BufWriter::new(file);
    let len = u64::try_from(index_ref.parameters.partition_number).unwrap(); // TODO expect
    tar_get::reserve_capacity(&mut writer, len).unwrap(); // TODO expect

    let file = writer.into_inner().unwrap();

    let file = Arc::new(Mutex::new(file));

    // for each partition, merge the corresponding bfs for every index
    for partition_idx in 0..index_ref.parameters.partition_number {
        // for each index, build the file name for the given partition
        let chunk_files: Vec<String> = indexes_metadata
            .iter()
            .enumerate()
            .map(|(index_dir, _)| {
                format!(
                    "{}/partition_bloom_filters_p{}.bin",
                    index_dir, partition_idx
                )
            })
            .collect();

        // collect the color counts from each index
        let color_counts: Vec<usize> = indexes_metadata.iter().map(|(_, count)| *count).collect();

        //merge
        merge_partition_bloom_filters(
            &chunk_files,
            partitioned_bf_size,
            index_ref.parameters.abundance_number.get(),
            &color_counts,
            &file,
        )?;
    }
    let parameters = Parameters {
        nb_color: new_color_number,
        ..index_ref.parameters
    };
    let index_merged = Reindeer2 {
        parameters,
        indexed_file_names,
        index_dir: String::from(output_dir),
    };
    index_merged.save_to_disk()?;

    log::info!(
        "Successfully merged {} indexes into directory: {}",
        indexes_metadata.len(),
        output_dir
    );
    Ok(())
}

// // --- MERGE FUNCTIONS ---

pub fn merge_all_partitions(
    chunk_files_dir: &str,
    output_dir: &str,
    partitioned_bf_size: usize,
    abundance_number: usize,
    color_counts_per_chunk: Vec<usize>, // number of colors in each chunk
    num_partitions: usize,
) -> io::Result<()> {
    let start_time = Instant::now();

    let path = Path::new(output_dir).join("partition_bloom_filters.bin");
    let file = OpenOptions::new()
        .read(true)
        .write(true)
        .create_new(true)
        .open(&path)
        .unwrap();
    let mut writer = BufWriter::new(file);
    let len = u64::try_from(num_partitions).expect("should have less than u64::MAX partitions");
    tar_get::reserve_capacity(&mut writer, len).expect("should have been able to write on disk");

    let file = writer.into_inner().unwrap();
    let file = Arc::new(Mutex::new(file));

    let nb_chunk = color_counts_per_chunk.len();
    // for each partition in parallel
    (0..num_partitions)
        // .into_par_iter() // TODO
        .try_for_each(|partition_idx| {
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
                &file,
            )?;

            Ok::<(), io::Error>(())
        })?;

    let elapsed_time = start_time.elapsed();
    log::info!(
        "All partitions merged and written to disk in {:.2?}",
        elapsed_time
    );

    Ok(())
}

pub fn merge_partition_bloom_filters(
    chunk_files: &[String],
    partitioned_bf_size: usize,
    abundance_number: usize,
    color_counts: &[usize], // number of colors for each chunk
    output_file: &Mutex<File>,
) -> io::Result<RoaringBitmap> {
    if chunk_files.is_empty() || chunk_files.len() != color_counts.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Mismatch between chunk files and color counts, or no chunks provided.",
        ));
    }

    // write all slices in the right order in a  larger filter
    let final_bf = merge_partition_slices_interleaved(
        chunk_files,
        partitioned_bf_size,
        abundance_number,
        color_counts,
    );

    let serializer =
        |writer: &mut BufWriter<&mut File>, bitmap: &RoaringBitmap| bitmap.serialize_into(writer);
    let mut output_file = output_file.lock().unwrap();
    tar_get::append_element(&mut output_file, &final_bf, serializer).unwrap();
    output_file.flush().unwrap();
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

fn merge_partition_slices_interleaved(
    chunk_files: &[String],
    partitioned_bf_size: usize,
    abundance_number: usize,
    color_counts: &[usize], // number of colors in each chunk
) -> RoaringBitmap {
    // preload all Bfs

    // OPTIMIZE we may be able to laod only one filter at a time
    let loaded_bfs: Vec<Option<(RoaringBitmap, usize)>> = chunk_files
        .iter()
        .map(|chunk_file| {
            match load_raw_bloom_filter(chunk_file) {
                Ok(bf) => Some(bf),
                Err(e) => {
                    log::error!("Failed to load Bloom filter {}: {}", chunk_file, e);
                    None // Return None for any errors
                }
            }
        })
        .collect();

    // vector to collect all positions for the final bf
    let mut final_positions = Vec::new();

    for slice_idx in 0..partitioned_bf_size {
        let mut current_offset = slice_idx * abundance_number * color_counts.iter().sum::<usize>();

        for (chunk_idx, loaded_bf) in loaded_bfs.iter().enumerate() {
            if let Some((chunk_bf, chunk_color_number)) = loaded_bf {
                let slice_start = slice_idx * abundance_number * chunk_color_number;
                let slice_end = slice_start + abundance_number * chunk_color_number;

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
                    chunk_bf
                        .range(slice_start_u32..slice_end_u32)
                        .map(|pos| pos - slice_start_u32 + current_offset_u32),
                );

                // update the offset for the next chunk in this slice
                current_offset += abundance_number * chunk_color_number;
            } else {
                log::error!("Skipping chunk {} due to previous load error", chunk_idx);
            }
        }
    }
    RoaringBitmap::from_sorted_iter(final_positions)
        .expect("Attempt to merge with unsorted positions")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_partition_bloom_filters() {
        use roaring::RoaringBitmap;
        use std::fs::{create_dir_all, File};
        use std::io::Write;

        let test_dir = "test_merge_partition_bloom_filters";
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

        let file = writer.into_inner().unwrap();

        let file = Arc::new(Mutex::new(file));

        // Call the modified function
        let merged_bf = merge_partition_bloom_filters(
            &chunk_files,
            partitioned_bf_size,
            abundance_number,
            &color_counts,
            &file,
        )
        .expect("Failed to merge partition Bloom filters");

        // check the merged Bloom filter
        // let merged_bf = merged_bloom_filter.lock().unwrap();

        //let final_output_path = format!("{}/partition_bloom_filters_p{}.bin", output_dir, partition_idx);
        //let (merged_bf, color_nb) = load_bloom_filter(&final_output_path).expect("Failed to load merged Bloom filter");

        //assert_eq!(
        //    color_nb, 5,
        //    "Expected the merged Bloom filter to represent 5 colors, but got {}",
        //    color_nb
        //);

        let expected_elements: Vec<u32> = vec![1, 2, 9, 10, 29];
        for elem in expected_elements {
            assert!(
                merged_bf.contains(elem),
                "Expected element {} not found in the merged Bloom filter",
                elem
            );
        }

        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");
    }
}
