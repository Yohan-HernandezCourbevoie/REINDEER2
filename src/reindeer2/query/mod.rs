mod format;

use std::cmp::min;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

#[derive(Copy, Clone, PartialEq, Debug)]
#[cfg_attr(any(test, debug_assertions), derive(Eq))] // ???
pub struct LogAbundance {
    value: u16,
}

impl LogAbundance {
    pub const QUERIED_BUT_ERROR: u16 = u16::MAX;
    pub const NOT_QUERIED: u16 = 0; // 0 so that hopefully compilers can optimize the creation of a slice of NOT_QUERIED
    pub const QUERIED_BUT_ABSENT: u16 = 1;
}

impl PartialOrd for LogAbundance {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for LogAbundance {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // removing 1 by wrapping
        // so that NOT_QUERIED is mapped to `u16::MAX`
        // so when we take the min, NOT_QUERIED is not preferred
        self.value.wrapping_sub(1).cmp(&other.value.wrapping_sub(1))
    }
}

impl LogAbundance {
    /// Accepts valid values from 1 to u8::MAX
    /// // 0 is mapped
    pub fn from_u8(value: u8) -> Self {
        Self {
            value: value as u16 + 1,
        }
    }

    /// Builds a `LogAbundance` from the result of a filter query.
    pub fn from_position_of_hit_in_the_filter(position_of_hit: u16) -> Self {
        // position of hit = 0 => present
        // as we have to encode "not queried" and "queried but absent", the first "present" must be 2
        let value = min(position_of_hit + 2, Self::QUERIED_BUT_ERROR - 1);
        Self { value }
    }

    pub fn new_error() -> Self {
        Self {
            value: Self::QUERIED_BUT_ABSENT,
        }
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct ApproxAbundance {
    value: u16,
}

impl PartialOrd for ApproxAbundance {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ApproxAbundance {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // removing 1 by wrapping
        // so that NOT_QUERIED is mapped to `u16::MAX`
        // so when we take the min, NOT_QUERIED is not preferred
        self.value.wrapping_sub(1).cmp(&other.value.wrapping_sub(1))
    }
}

impl ApproxAbundance {
    pub const QUERIED_BUT_ERROR: u16 = LogAbundance::QUERIED_BUT_ERROR;
    pub const NOT_QUERIED: u16 = LogAbundance::NOT_QUERIED;
    pub const QUERIED_BUT_ABSENT: u16 = LogAbundance::QUERIED_BUT_ABSENT;

    // pub fn from_dense(dense_result: u8, base: f64) -> Self {
    //     let dense_result = dense_result as u16 + 1;
    //     let approx_count = if dense_result == LogAbundance::QUERIED_BUT_ERROR {
    //         Self::QUERIED_BUT_ERROR
    //     } else if dense_result == LogAbundance::NOT_QUERIED {
    //         Self::NOT_QUERIED
    //     } else if dense_result == LogAbundance::QUERIED_BUT_ABSENT {
    //         Self::QUERIED_BUT_ABSENT
    //     } else {
    //         dbg!(approximate_value(dense_result - 2, base) + 2) // FIXME limit this
    //     };
    //     Self {
    //         value: approx_count,
    //     }
    // }

    pub fn from_log(log_abund: LogAbundance, base: f64) -> Self {
        let approx_count = if log_abund.value == LogAbundance::QUERIED_BUT_ERROR {
            Self::QUERIED_BUT_ERROR
        } else if log_abund.value == LogAbundance::NOT_QUERIED {
            Self::NOT_QUERIED
        } else if log_abund.value == LogAbundance::QUERIED_BUT_ABSENT {
            Self::QUERIED_BUT_ABSENT
        } else {
            approximate_value(log_abund.value - 2, base) + 2 // FIXME limit this
        };
        Self {
            value: approx_count,
        }
    }

    pub fn from_position_of_hit_in_the_filter(hit_position: u16, base: f64) -> Self {
        let value = approximate_value(hit_position, base) + 2; // FIXME limit this
        Self { value }
    }

    pub fn new(val: u16) -> Self {
        Self { value: val + 2 }
    }

    pub fn new_not_queried() -> Self {
        Self {
            value: Self::NOT_QUERIED,
        }
    }

    const fn is_valid_abundance(&self) -> bool {
        (self.value != Self::QUERIED_BUT_ERROR) && (self.value != Self::NOT_QUERIED)
    }

    // const fn is_valid_and_absent(&self) -> bool {
    //     if self.is_valid_abundance() {
    //         self.value == LogAbundance::QUERIED_BUT_ABSENT
    //     } else {
    //         false
    //     }
    // }

    const fn is_not_queried(&self) -> bool {
        self.value == LogAbundance::NOT_QUERIED
    }

    pub const fn to_value(&self) -> Option<u16> {
        if self.is_valid_abundance() {
            if self.value == LogAbundance::QUERIED_BUT_ABSENT {
                Some(0)
            } else {
                Some(self.value - 2) // FIXME +1 ? -1 ?
            }
        } else {
            None
        }
    }

    pub fn select_abundance_from_candidates(candidates: &[(u32, Self)]) -> Option<&(u32, Self)> {
        candidates.iter().min()
    }

    pub fn new_error() -> Self {
        Self {
            value: Self::QUERIED_BUT_ERROR,
        }
    }
}

use super::{
    approximate_value, compute_base_position,
    dense_index::DenseIndexPartition,
    load_bloom_filter,
    minimizer_iter::{kmer_minimizers_sampled, Sampler},
};
use bio::io::fasta;
pub use format::{write_header, write_kmer_query, EnrichedOutputFormat};

/// Builds a map from partition_index -> Vec of (sequence_id, position_kmer_in_sequence, kmer_hash).
pub fn build_partitions_kmers<S>(
    batch: &[fasta::Record],
    k: usize,
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
        let kmer_minimizers = kmer_minimizers_sampled(sequence, k, m, canonical, sampler)
            .expect("should have been able to iterate over kmers");

        for (position, (kmer_hash, minimizer)) in kmer_minimizers.enumerate() {
            let partition_index = (minimizer % partition_number) as usize;
            partition_kmers.entry(partition_index).or_default().push((
                record_id,
                position as u32,
                kmer_hash,
            ));
        }
    }

    partition_kmers
}

// TODO rename
/// parameter `kmers` is Vec<(read_id, pos_in_read, hash)>
pub fn fold_into_hashmap(
    mut local_results: HashMap<usize, Vec<Vec<(u32, ApproxAbundance)>>>,
    partition_index: usize,
    kmers: Vec<(usize, u32, u64)>,
    base: f64,
    bf_dir: &str,
    bf_size: u64,
    partition_number: usize,
    color_number: usize,
    abundance_number: usize,
    is_dense: bool,
) -> HashMap<usize, Vec<Vec<(u32, ApproxAbundance)>>> {
    // Load the partition's Bloom filter
    let path_bf = format!(
        "{}/partition_bloom_filters_p{}.bin",
        bf_dir, partition_index
    );
    let maybe_bf = load_bloom_filter(&path_bf);

    if let Ok((bitmap, _maybe_aux_data)) = maybe_bf {
        let hashmap: DenseIndexPartition = if is_dense {
            let path_dense_index =
                format!("{}/partition_dense_index_p{}.bin", bf_dir, partition_index);
            DenseIndexPartition::load_from_disk(&path_dense_index).expect(&format!(
                "Failed to load dense index for partition {}",
                partition_index
            ))
        } else {
            DenseIndexPartition::new()
        };
        //  For each k-mer in this partition
        for (sequence_id, kmer_position, kmer_hash) in kmers {
            let color_abundances = if hashmap.contains_key(&kmer_hash) {
                let log_abundance_vector = hashmap
                    .get_abundance(&kmer_hash)
                    .expect("failed to read the hashmap");
                // OPTIMIZE use only a vector instead of a vec of vec
                let mut color_abundances = vec![Vec::new(); color_number];
                for (color, log_abundance) in log_abundance_vector.iter().enumerate() {
                    color_abundances[color]
                        .push((kmer_position, LogAbundance::from_u8(*log_abundance)));
                }
                color_abundances
            } else {
                // Compute base position
                let base_position = compute_base_position(
                    kmer_hash,
                    (bf_size as usize) / partition_number,
                    color_number,
                    abundance_number,
                );

                // color_abundances[color] -> Vec of (log) counts for that color
                let mut color_abundances = vec![Vec::new(); color_number];
                // TODO this function could have a better name
                update_color_abundances(
                    &bitmap,
                    base_position,
                    color_number,
                    abundance_number,
                    kmer_position,
                    base,
                    &mut color_abundances,
                );
                color_abundances
            };

            // Convert log abundances to approximate integer counts
            let approximate_counts: Vec<Vec<(u32, ApproxAbundance)>> = color_abundances
                .into_iter()
                .map(|abunds_for_color| {
                    abunds_for_color
                        .into_iter()
                        .map(|(kmer_pos, log_abund)| {
                            (kmer_pos, ApproxAbundance::from_log(log_abund, base))
                        })
                        .collect()
                })
                .collect();

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

// TOUN
// TODO why the outparameter ?
/// update color abundances for a specific base position in the Bloom filter
pub fn update_color_abundances(
    bitmap: &roaring::RoaringBitmap,
    base_position: u64,
    color_number: usize,
    abundance_number: usize,
    kmer_position: u32,
    base: f64,
    color_abundances: &mut [Vec<(u32, LogAbundance)>],
) {
    debug_assert!(abundance_number < u16::MAX as usize); // TODO make this check before ?
    let abundance_number = abundance_number as u16; // TODO take u16 as parameter ?
    for color in 0..color_number {
        let mut insert = false;
        for abundance in 0..abundance_number {
            let position_to_check =
                base_position + (color as u64) * (abundance_number as u64) + (abundance as u64);

            if bitmap.contains(position_to_check as u32) {
                color_abundances[color].push((
                    kmer_position,
                    LogAbundance::from_position_of_hit_in_the_filter(abundance),
                ));
                insert = true;
                break; // keep the minimum
            }
        }
        if !insert {
            // TODO weird discuss value
            color_abundances[color].push((kmer_position, LogAbundance::new_error()));
            // important to record absent k-mers, to compute the median value, also, todo test
        }
    }
}

fn load_kmer_counts_vector(dir_path: &str) -> io::Result<Vec<usize>> {
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
    abund_values: Vec<(u32, ApproxAbundance)>,
    nb_kmer_in_query: usize,
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
            debug_assert!((**max as usize) < set_positions.len());
        }
        debug_assert_eq!(
            **set_positions.iter().max().unwrap() as usize,
            abund_values.len() - 1
        );
    }

    let mut abund_values_ordered = vec![ApproxAbundance::new_not_queried(); nb_kmer_in_query];
    for (kmer_pos, abundance) in abund_values {
        abund_values_ordered[kmer_pos as usize] = abundance;
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
    fn approx_back_and_forth() {
        for i in [0, 4, 8, 9, 4, 456, 789, 54] {
            let approx = ApproxAbundance::new(i);
            assert_eq!(approx.to_value().unwrap(), i);
        }
    }
    #[test]
    fn test_update_color_abundances() {
        use roaring::RoaringBitmap;

        let mut bitmap = RoaringBitmap::new();
        let base_position = 100;
        let color_number = 3;
        let abundance_number = 2;
        let base = 2.0;

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
            vec![(0, LogAbundance::from_u8(1))], // color 0 has abundance levels 0 and 1 -> will keep the min
            vec![(0, LogAbundance::from_u8(1))], // color 1 has abundance level 0
            vec![(0, LogAbundance::new_error())], // color 2 has no abundance, set to QUERIED_BUT_ERROR instead (handle later in the pipeline)
        ];

        assert_eq!(
            color_abundances, expected_color_abundances,
            "Color abundances mismatch: expected {:?}, got {:?}",
            expected_color_abundances, color_abundances
        );
    }
}
