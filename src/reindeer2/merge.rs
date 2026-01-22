use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::{
    fs::{self, File},
    io::{self, BufRead, BufWriter, Write},
    time::Instant,
};

use crate::reindeer2::{load_bloom_filter, Parameters, Reindeer2};

// merge an arbitrary number of indexes  listed in a fof
// each line of the fof is expected to he path to one index dir
// in the end, alsos write a  new metadata JSON file in the output dir
pub fn merge_multiple_indexes(indexes_fof: &str, output_dir: &str) -> io::Result<()> {
    // read the list of index directories.
    let index_dirs: Vec<String> = {
        let file = File::open(indexes_fof)?;
        io::BufReader::new(file)
            .lines()
            .filter_map(Result::ok)
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
        let mut merged_bf = RoaringBitmap::new();
        merge_partition_bloom_filters(
            chunk_files,
            partition_idx,
            partitioned_bf_size,
            index_ref.parameters.abundance_number.get(),
            &color_counts,
            &mut merged_bf,
            output_dir,
            new_color_number,
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
    total_nb_colors: usize,
) -> io::Result<()> {
    let start_time = Instant::now();
    // for each partition in parallel
    (0..num_partitions)
        .into_par_iter()
        .try_for_each(|partition_idx| {
            // collect chunk files for the current partition
            let chunk_files_for_partition: Vec<String> = color_counts_per_chunk
                .iter()
                .enumerate()
                .map(|(chunk_idx, _)| {
                    format!(
                        "{}/partition_bloom_filters_c{}_p{}.bin",
                        chunk_files_dir, chunk_idx, partition_idx
                    )
                })
                .collect();

            // bf for this partition
            let mut partition_bf = RoaringBitmap::new();

            // merge all chunks for the current partition + serialize
            merge_partition_bloom_filters(
                chunk_files_for_partition,
                partition_idx,
                partitioned_bf_size,
                abundance_number,
                &color_counts_per_chunk,
                &mut partition_bf,
                output_dir,
                total_nb_colors,
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
    chunk_files: Vec<String>,
    partition_idx: usize,
    partitioned_bf_size: usize,
    abundance_number: usize,
    color_counts: &[usize], // number of colors for each chunk
    bloom_filter: &mut RoaringBitmap,
    output_dir: &str,
    total_nb_colors: usize,
) -> io::Result<()> {
    if chunk_files.is_empty() || chunk_files.len() != color_counts.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Mismatch between chunk files and color counts, or no chunks provided.",
        ));
    }

    // write all slices in the right order in a  larger filter
    let final_bf = merge_partition_slices_interleaved(
        &chunk_files,
        partitioned_bf_size,
        abundance_number,
        color_counts,
    );
    let output_file_path = format!(
        "{}/partition_bloom_filters_p{}.bin",
        output_dir, partition_idx
    );
    let output_file = File::create(&output_file_path)?;
    let mut writer = BufWriter::new(output_file);
    writer.write_all(&total_nb_colors.to_le_bytes())?;
    final_bf.serialize_into(&mut writer)?;
    // update the bloom filter for the partition
    *bloom_filter = final_bf;

    Ok(())
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
    let loaded_bfs: Vec<Option<(RoaringBitmap, usize)>> = chunk_files
        .iter()
        .map(|chunk_file| {
            match load_bloom_filter(chunk_file) {
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
