mod format;

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Read};
use std::path::Path;

use super::{
    approximate_value, compute_base, compute_base_position, dense_index::DenseIndexPartition,
    kmer_minimizers_seq_level, load_bloom_filter, process_fasta_in_batches, read_file,
    read_indexed_file_names, update_color_abundances,
};
use crate::reindeer2::OutputFormat;
use bio::io::fasta;
use format::{write_header, write_kmer_query, EnrichedOutputFormat};
use rayon::prelude::*;

// === QUERY ===

// --- MAIN FUNCTION ---

// --- QUERY FUNCTIONS ---

/// Builds a map from partition_index -> Vec of (sequence_id, position_kmer_in_sequence, kmer_hash).
fn build_partitions_kmers(
    batch: &[fasta::Record],
    k: usize,
    m: usize,
    partition_number: u64,
    canonical: bool,
) -> HashMap<usize, Vec<(usize, usize, u64)>> {
    let mut partition_kmers: HashMap<usize, Vec<(usize, usize, u64)>> = HashMap::new();
    for (record_id, record) in batch.iter().enumerate() {
        let kmer_minimizers = kmer_minimizers_seq_level(record.seq(), k, m, canonical);
        for (position, (kmer_hash, minimizer)) in kmer_minimizers.enumerate() {
            let partition_index = (minimizer % partition_number) as usize;
            partition_kmers
                .entry(partition_index)
                .or_default()
                .push((record_id, position, kmer_hash));
        }
    }
    partition_kmers
}

// TODO rename
/// parameter `kmers` is Vec<(read_id, pos_in_read, hash)>
fn fold_into_hashmap(
    mut local_results: HashMap<usize, Vec<Vec<(usize, u16)>>>,
    partition_index: usize,
    kmers: Vec<(usize, usize, u64)>,
    base: f64,
    bf_dir: &str,
    bf_size: u64,
    partition_number: usize,
    color_number: usize,
    abundance_number: usize,
    is_dense: bool,
) -> HashMap<usize, Vec<Vec<(usize, u16)>>> {
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
                let mut color_abundances = vec![Vec::new(); color_number];
                for (color, log_abundance) in log_abundance_vector.iter().enumerate() {
                    if *log_abundance == 0 {
                        // TODO discuss weird value
                        color_abundances[color].push((kmer_position, 666));
                    } else {
                        color_abundances[color]
                            .push((kmer_position, (*log_abundance - 1) as usize));
                    }
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
                    &mut color_abundances,
                );
                color_abundances
            };

            // Convert log abundances to approximate integer counts
            let approximate_counts: Vec<Vec<(usize, u16)>> = color_abundances
                .into_iter()
                .map(|abunds_for_color| {
                    abunds_for_color
                        .into_iter()
                        .map(|(kmer_pos, log_abund)| {
                            let approx_count = if log_abund == 666 {
                                0
                            } else {
                                approximate_value(log_abund, base)
                            };
                            (kmer_pos, approx_count)
                        })
                        .collect()
                })
                .collect();

            // Accumulate results in local_results
            let entry = local_results
                .entry(sequence_id)
                .or_insert_with(|| vec![Vec::new(); color_number]);
            for (color_idx, approx_values) in approximate_counts.into_iter().enumerate() {
                entry[color_idx].push(
                    *approx_values
                        .iter()
                        .min()
                        .expect("An abundance vector returned empty"),
                );
            }
        }
    } else {
        eprintln!(
            "Failed to load Bloom filter for partition {}",
            partition_index
        );
    }

    local_results
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

/// Formats a fasta header by removing the first `>` and taking up to the first space (excluded).
/// E.g.: `>seq1 ka:f:30` -> `seq1`
fn strip_header(s: &str) -> &str {
    let stripped = if let Some(stripped) = s.strip_prefix('>') {
        stripped
    } else {
        s
    };
    stripped.split(' ').next().unwrap()
}

fn sort_abundance_vec(abund_values: Vec<(usize, u16)>) -> Vec<u16> {
    use itertools::Itertools;

    #[cfg(debug_assertions)]
    {
        use std::collections::HashSet;

        let set_positions: HashSet<&usize> = abund_values
            .iter()
            .map(|(kmer_pos, _abundance)| kmer_pos)
            .collect();
        debug_assert_eq!(set_positions.len(), abund_values.len());
        debug_assert_eq!(
            **set_positions.iter().max().unwrap(),
            abund_values.len() - 1
        );
        debug_assert_eq!(**set_positions.iter().min().unwrap(), 0);
    }
    abund_values
        .into_iter()
        .sorted_by_key(|(kmer_pos, _abundance)| *kmer_pos)
        .map(|(_kmer_pos, abundance)| abundance)
        .collect()
}

// TODO use compile time to decide wether to sort the output or not
fn query_single_fasta_batch(
    bf_dir: &str,
    k: usize,
    m: usize,
    bf_size: u64,
    partition_number: usize,
    color_number: usize,
    abundance_number: usize,
    dense_option: bool,
    canonical: bool,
    base: f64,
    batch: &[fasta::Record],
) -> HashMap<usize, Vec<Vec<u16>>> {
    let partition_kmers = build_partitions_kmers(batch, k, m, partition_number as u64, canonical);

    // --- PARALLEL PHASE: process each partition's k-mers in parallel ---
    // This will store final results for *all* sequences in this batch.
    // Key: sequence header; Value: vector of vector of abundances
    let result_with_positions: HashMap<usize, Vec<Vec<(usize, u16)>>> = partition_kmers
        .into_par_iter()
        // 1) Create a local HashMap in each thread
        .fold(
            // OPTIMIZE use compile time monomorphization to get rid of this (usize, u16)
            // (I expect to go from 128 per element in the vector to 16)
            // (I needed this extra usize to get the correct order of kmer in the output)
            HashMap::<usize, Vec<Vec<(usize, u16)>>>::new,
            |local_results: HashMap<usize, Vec<Vec<(usize, u16)>>>, (partition_index, kmers)| {
                fold_into_hashmap(
                    local_results,
                    partition_index,
                    kmers,
                    base,
                    bf_dir,
                    bf_size,
                    partition_number,
                    color_number,
                    abundance_number,
                    dense_option,
                )
            },
        )
        // 2) Reduce all local HashMaps into a single HashMap
        .reduce(HashMap::<usize, Vec<Vec<(usize, u16)>>>::new, merge_results);

    result_with_positions
        .into_iter()
        .map(|(header, color_vectors)| {
            (
                header,
                color_vectors.into_iter().map(sort_abundance_vec).collect(),
            )
        })
        .collect()
}

pub fn query_sequences_in_batches(
    fasta_file: &str,
    bf_dir: &str,
    k: usize,
    m: usize,
    bf_size: u64,
    partition_number: usize,
    color_number: usize,
    abundance_number: usize,
    abundance_max: u16,
    batch_size: usize,
    output_file: &str,
    dense_option: bool,
    output_format: OutputFormat,
    coverage: f32,
    canonical: bool,
) -> io::Result<()> {
    let reader = read_file(fasta_file)?;
    let mut writer = BufWriter::new(File::create(output_file)?);
    let output_format = EnrichedOutputFormat::from_pub_output_format(output_format, bf_dir);
    // write the header of the result file
    write_header(bf_dir, &output_format, &mut writer)
        .expect("should have been able to write the header of the result file");

    let base = compute_base(abundance_number, abundance_max);
    let output_format = &output_format;

    // Process FASTA in chunks of `batch_size` records
    process_fasta_in_batches(reader, batch_size, |batch| {
        let sequence_results = query_single_fasta_batch(
            bf_dir,
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            dense_option,
            canonical,
            base,
            batch,
        );
        // Now `sequence_results` has the combined data for this batch.
        // Let's compute the output in the requested format.
        write_kmer_query(
            batch,
            output_format,
            coverage,
            sequence_results,
            &mut writer,
        )
        .expect("should have been able to write the query result");
    })
    .expect("should have been able to process fasta files");

    Ok(())
}

fn merge_results(
    mut acc: HashMap<usize, Vec<Vec<(usize, u16)>>>,
    local: HashMap<usize, Vec<Vec<(usize, u16)>>>,
) -> HashMap<usize, Vec<Vec<(usize, u16)>>> {
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
