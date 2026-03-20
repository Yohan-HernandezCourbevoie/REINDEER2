use itertools::Itertools;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::{
    cmp::min,
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
    time::Instant,
};

use crate::reindeer2::{
    create_and_reserve_tar_get_file, filter::load_raw_bloom_filter,
    merge::merge_partition_slices_interleaved, NB_FILE_IN_AN_INDEX,
};

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
        .into_par_iter()
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
fn merge_partition_bloom_filters(
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
