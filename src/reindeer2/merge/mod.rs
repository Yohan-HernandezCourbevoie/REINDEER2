mod merge_indexes;
mod merge_partitions;

use roaring::RoaringBitmap;

pub use merge_indexes::merge_multiple_indexes;
pub use merge_partitions::{merge_all_partitions_of_chunks, remove_merged_partitions};

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
