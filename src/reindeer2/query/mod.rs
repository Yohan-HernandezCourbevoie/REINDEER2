mod approx_abundance;
mod findere;
mod format;

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

use crate::reindeer2::NB_FILE_IN_AN_INDEX;

use super::{
    compute_base_position,
    dense_index::DenseIndexPartition,
    filter::load_bloom_filter_from_big_file,
    minimizer_iter::{kmer_minimizers_all, Sampler},
};
use bio::io::fasta;
pub use format::{write_header, write_kmer_query, EnrichedOutputFormat};

pub use approx_abundance::ApproxAbundance;
pub use findere::fimpera;

/// Builds a map from partition_index -> Vec of (sequence_id, position_kmer_in_sequence, smer_hash).
pub fn build_partitions_smers<S>(
    batch: &[fasta::Record],
    s: usize,
    m: usize,
    partition_number: u64,
    canonical: bool,
    sampler: &S,
) -> HashMap<usize, Vec<(usize, u32, u64)>>
where
    S: Sampler,
{
    let mut partition_kmers: HashMap<usize, Vec<(usize, u32, u64)>> = HashMap::new();
    for (record_id, record) in batch.iter().enumerate() {
        let sequence = record.seq();
        // covers the case "usize is more than 64 bits"
        assert!(
            sequence.len() < (u32::MAX as usize),
            "queried sequence length should be smaller than 2^32"
        );
        // covers the case "usize is less than 64 bits"
        assert!(
            (sequence.len() as u32) < u32::MAX,
            "queried sequence length should be smaller than 2^32"
        );
        let kmer_minimizers = kmer_minimizers_all(sequence, s, m, canonical)
            .expect("should have been able to iterate over kmers");

        for (position, (smer_hash, minimizer)) in kmer_minimizers.enumerate() {
            if sampler.filter((smer_hash, minimizer)) {
                let partition_index = (minimizer % partition_number) as usize;
                partition_kmers.entry(partition_index).or_default().push((
                    record_id,
                    position as u32,
                    smer_hash,
                ));
            }
        }
    }

    partition_kmers
}

// TODO rename
/// parameter `smers` is Vec<(read_id, pos_in_read, hash)>
pub fn fold_into_hashmap(
    mut local_results: HashMap<usize, Vec<Vec<(u32, ApproxAbundance)>>>,
    partition_index: usize,
    smers: Vec<(usize, u32, u64)>,
    base: f64,
    bf_dir: &str,
    bf_size: u64,
    partition_number: usize,
    color_number: usize,
    abundance_number: usize,
    is_dense: bool,
) -> HashMap<usize, Vec<Vec<(u32, ApproxAbundance)>>> {
    // Load the partition's Bloom filter
    let nb_partition_in_a_file = partition_number.div_ceil(NB_FILE_IN_AN_INDEX);
    let group = partition_index / nb_partition_in_a_file;
    let index = partition_index % nb_partition_in_a_file;
    let path_bf = format!("{}/partition_bloom_filters_group{}.bin", bf_dir, group);
    let maybe_bf = load_bloom_filter_from_big_file(&path_bf, index as u64); // TODO conversion error, expect

    if let Ok(bitmap) = maybe_bf {
        let hashmap: DenseIndexPartition = if is_dense {
            let path_dense_index =
                format!("{}/partition_dense_index_p{}.bin", bf_dir, partition_index);
            DenseIndexPartition::load_from_disk(&path_dense_index).unwrap_or_else(|_| {
                panic!(
                    "Failed to load dense index for partition {}",
                    partition_index
                )
            })
        } else {
            DenseIndexPartition::new()
        };
        //  For each k-mer in this partition
        for (sequence_id, smer_position, smer_hash) in smers {
            let approximate_counts = if hashmap.contains_key(&smer_hash) {
                let log_abundance_vector = hashmap
                    .get_abundance(&smer_hash)
                    .expect("failed to read the hashmap");
                // OPTIMIZE use only a vector instead of a vec of vec
                let mut approximate_counts = vec![Vec::new(); color_number];
                for (color, log_abundance) in log_abundance_vector.iter().enumerate() {
                    approximate_counts[color].push((
                        smer_position,
                        ApproxAbundance::from_dense(*log_abundance, base),
                    ));
                }
                approximate_counts
            } else {
                // Compute base position
                let base_position = compute_base_position(
                    smer_hash,
                    (bf_size as usize) / partition_number,
                    color_number,
                    abundance_number,
                );

                // color_abundances[color] -> Vec of (log) counts for that color
                let mut approximate_counts = vec![Vec::new(); color_number];
                // TODO this function could have a better name
                update_color_abundances(
                    &bitmap,
                    base_position,
                    color_number,
                    abundance_number,
                    smer_position,
                    base,
                    &mut approximate_counts,
                );
                approximate_counts
            };

            // Accumulate results in local_results
            let entry = local_results
                .entry(sequence_id)
                .or_insert_with(|| vec![Vec::new(); color_number]);
            // FIXME: remove, maybe replace by something else ?
            for (color_idx, approx_values) in approximate_counts.into_iter().enumerate() {
                entry[color_idx].push(
                    // FIXME: don't we take one the smallest first element ?
                    *ApproxAbundance::select_abundance_from_candidates(&approx_values)
                        .expect("An abundance vector returned empty"),
                );
            }
        }
    } else {
        log::error!(
            "Failed to load Bloom filter for partition {}",
            partition_index
        );
    }

    local_results
}

// TODO this should be in the Filter
fn query_smer(
    bitmap: &roaring::RoaringBitmap,
    base_position: u64,
    abundance_number: u16,
    base: f64,
    color: usize,
) -> ApproxAbundance {
    let mut res = ApproxAbundance::new_absent();
    for abundance in 0..abundance_number {
        let position_to_check =
            base_position + (color as u64) * (abundance_number as u64) + (abundance as u64);

        if bitmap.contains(position_to_check as u32) {
            res = ApproxAbundance::from_position_of_hit_in_the_filter(abundance, base);
        }
    }
    res
}

// TOUN
// TODO why the outparameter ?
/// update color abundances for a specific base position in the Bloom filter
pub fn update_color_abundances(
    bitmap: &roaring::RoaringBitmap,
    base_position: u64,
    color_number: usize,
    abundance_number: usize,
    smer_position: u32,
    base: f64,
    color_abundances: &mut [Vec<(u32, ApproxAbundance)>],
) {
    debug_assert!(abundance_number < u16::MAX as usize); // TODO make this check before ?
    let abundance_number: u16 = abundance_number as u16; // TODO take u16 as parameter ?
    for color in 0..color_number {
        let smer_approx_abundance =
            query_smer(bitmap, base_position, abundance_number, base, color);
        color_abundances[color].push((smer_position, smer_approx_abundance));
    }
}

pub fn load_kmer_counts_vector(dir_path: &str) -> io::Result<Vec<usize>> {
    let mut file = File::open(Path::new(dir_path).join("kmer_counts_per_color.bin"))?;

    // Read the rest of the file to deserialize the hashmap
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;
    let counts_vector = bincode::deserialize_from(&buffer[..]).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "Failed to deserialize the counts vector",
        )
    })?;
    Ok(counts_vector)
}

// /// Formats a fasta header by removing the first `>` and taking up to the first space (excluded).
// /// E.g.: `>seq1 ka:f:30` -> `seq1`
// fn strip_header(s: &str) -> &str {
//     let stripped = if let Some(stripped) = s.strip_prefix('>') {
//         stripped
//     } else {
//         s
//     };
//     stripped.split(' ').next().unwrap()
// }

pub fn sort_abundance_vec(
    abund_values: &Vec<(u32, ApproxAbundance)>,
    nb_smer_in_query: usize,
    #[cfg(any(debug_assertions, test))] s: usize,
) -> Vec<ApproxAbundance> {
    #[cfg(debug_assertions)]
    {
        use std::collections::HashSet;

        let set_positions: HashSet<&u32> = abund_values
            .iter()
            .map(|(kmer_pos, _abundance)| kmer_pos)
            .collect();
        // there are no duplicates in the positions of abund_values
        debug_assert_eq!(set_positions.len(), abund_values.len());
        // all positions are in the range of possible positions
        // FIXME what if usize is less than 32 bits ?
        if let Some(max) = set_positions.iter().max() {
            debug_assert!((**max as usize) < nb_smer_in_query + s - 1);
        }
    }

    let mut abund_values_ordered = vec![ApproxAbundance::new_not_queried(); nb_smer_in_query];
    for (kmer_pos, abundance) in abund_values {
        abund_values_ordered[*kmer_pos as usize] = *abundance;
    }

    abund_values_ordered
}

pub fn merge_results(
    mut acc: HashMap<usize, Vec<Vec<(u32, ApproxAbundance)>>>,
    local: HashMap<usize, Vec<Vec<(u32, ApproxAbundance)>>>,
) -> HashMap<usize, Vec<Vec<(u32, ApproxAbundance)>>> {
    for (seq_id, color_vecs) in local {
        let entry = acc
            .entry(seq_id)
            .or_insert_with(|| vec![Vec::new(); color_vecs.len()]);

        // Merge each color's abundances
        for (color_idx, local_abunds) in color_vecs.iter().enumerate() {
            entry[color_idx].extend(local_abunds.iter().copied());
        }
    }
    acc
}

// TODO write more unit tests for this file
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_update_color_abundances() {
        use roaring::RoaringBitmap;

        let mut bitmap = RoaringBitmap::new();
        let base_position = 100;
        let color_number = 3;
        let abundance_number = 2;
        let base = 2.0;

        #[allow(clippy::identity_op)]
        bitmap.insert((base_position + 0) as u32); // color 0, abundance 0
        bitmap.insert((base_position + 1) as u32); // color 0, abundance 1
        bitmap.insert((base_position + 2) as u32); // color 1, abundance 0

        let mut color_abundances = vec![vec![]; color_number];

        update_color_abundances(
            &bitmap,
            base_position,
            color_number,
            abundance_number,
            0,
            base,
            &mut color_abundances,
        );

        let expected_color_abundances = vec![
            vec![(0, ApproxAbundance::from_dense(2, base))], // color 0 has abundance levels 0 and 1 -> will keep the max
            vec![(0, ApproxAbundance::from_dense(1, base))], // color 1 has abundance level 0
            vec![(0, ApproxAbundance::new_absent())], // color 2 has no abundance, set to QUERIED_BUT_ERROR instead (handle later in the pipeline)
        ];

        assert_eq!(
            color_abundances, expected_color_abundances,
            "Color abundances mismatch: expected {:?}, got {:?}",
            expected_color_abundances, color_abundances
        );
    }
}
