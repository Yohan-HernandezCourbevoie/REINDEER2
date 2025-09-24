use bio::io::fasta;
use roaring::RoaringBitmap;
use std::collections::HashMap;
use std::io;
use std::sync::{atomic, Arc, Mutex};

use crate::reindeer2::{
    compute_location_filter, compute_log_abundance, extract_count, kmer_minimizers_seq_level,
    process_fasta_in_batches, read_file, HeaderType,
};

// --- INDEX FUNCTIONS ---

pub fn process_fasta_file(
    path: &str,
    maybe_dense_indexes: &Option<Arc<Vec<Mutex<HashMap<u64, Vec<u8>>>>>>,
    bloom_filters: &Arc<Vec<Mutex<RoaringBitmap>>>,
    k: usize,
    m: usize,
    partitioned_bf_size: usize,
    partition_number: usize,
    color_number: usize,
    color_number_global: usize,
    threshold: usize,
    abundance_number: usize,
    abundance_max: u16,
    path_num: usize,
    path_num_global: usize,
    base: f64,
    chunk_index: usize,
    header_type: HeaderType,
    max_map_size: usize, // Maximum size for the hash map
    total_kmers: &atomic::AtomicU64,
    atomic_dense_kmers_count: &atomic::AtomicU64,
    atomic_sparse_kmers_count: &atomic::AtomicU64,
    kmer_counts_vector: &Arc<Mutex<Vec<usize>>>,
    canonical: bool,
) -> io::Result<()> {
    let atomic_record_count = atomic::AtomicU64::new(0);
    let mut kmer_count: usize = 0;
    let reader = read_file(path)?;

    // this part fills the BFs per partition
    fn flush_map(
        partition_kmers: &mut HashMap<usize, Vec<(u64, u16, usize, usize)>>,
        partitioned_bf_size: usize,
        color_number: usize,
        abundance_number: usize,
        bloom_filters: &Arc<Vec<Mutex<RoaringBitmap>>>,
    ) {
        for (partition_index, kmers) in partition_kmers.drain() {
            // iterates and empties the hash map when needed

            let mut kmer_hashes_to_update = Vec::with_capacity(kmers.len());
            for (kmer_hash, log_abundance, path_num, _chunk_index) in kmers {
                // select the bit in the BF
                let position = compute_location_filter(
                    kmer_hash,
                    partitioned_bf_size,
                    color_number,
                    path_num,
                    abundance_number,
                    log_abundance,
                );
                kmer_hashes_to_update.push(position as u32); // accumulate bits to be modified for this bf
            }
            let mut bloom_filter = bloom_filters[partition_index] // select the correct BF for the given partition
                .lock()
                .expect("Failed to lock bloom filter");
            // TODO this takes a lot of time
            bloom_filter.extend(kmer_hashes_to_update);
        }
    }

    process_fasta_in_batches(reader, 10_000, |batch| {
        // read file 10_000 lines at once
        let mut partition_kmers: HashMap<usize, Vec<(u64, u16, usize, usize)>> = HashMap::new(); // keep kmer info to fill BFs

        // this part reads the batches of sequences and records kmers info until the structure is too large in memory
        for record in batch {
            let processed = process_fasta_record(Ok(record), base, abundance_max, &header_type); //read fasta
            match processed {
                Ok((seq, log_abundance, count_value)) => {
                    if log_abundance != 666 {
                        // case where the abundance value of the kmers in the unitigs file was < 1
                        atomic_record_count.fetch_add(1, atomic::Ordering::Relaxed);
                        let seq_str = std::str::from_utf8(&seq).expect("Invalid UTF-8 sequence");

                        for (kmer_hash, minimizer) in
                            kmer_minimizers_seq_level(seq_str.as_bytes(), k, m, canonical)
                        {
                            kmer_count += count_value as usize;
                            //for (kmer_hash, (minimizer, _)) in nt_hash_iterator.zip(min_iter) { // iterate on both minimizer and hash for each kmer
                            let partition_index = (minimizer % (partition_number as u64)) as usize;
                            total_kmers.fetch_add(1, atomic::Ordering::Relaxed);

                            match maybe_dense_indexes {
                                Some(dense_indexes) => {
                                    let mut dense_index = dense_indexes[partition_index]
                                        .lock()
                                        .expect("Failed to lock the dense index");

                                    // write in the dense index if the k-mer can be dense, else, put it in the hashmap for sparses
                                    if dense_index.contains_key(&kmer_hash) {
                                        // update the vector with the right abundance
                                        if let Some(abundance_vector) =
                                            dense_index.get_mut(&kmer_hash)
                                        {
                                            abundance_vector[path_num_global] =
                                                (log_abundance + 1) as u8;
                                        }
                                        atomic_dense_kmers_count
                                            .fetch_add(1, atomic::Ordering::Relaxed);
                                    } else if path_num_global <= threshold {
                                        // create a new abundance vector for the k-mer
                                        let mut abundance_vector: Vec<u8> =
                                            vec![0; color_number_global];
                                        abundance_vector[path_num_global] =
                                            (log_abundance + 1) as u8;
                                        dense_index.insert(kmer_hash, abundance_vector);
                                        atomic_dense_kmers_count
                                            .fetch_add(1, atomic::Ordering::Relaxed);
                                    } else {
                                        // write the k-mer in a file of sparse k-mer from this color
                                        partition_kmers // separate the kmers per partition
                                            .entry(partition_index)
                                            .or_default()
                                            .push((
                                                kmer_hash,
                                                log_abundance,
                                                path_num,
                                                chunk_index,
                                            ));
                                        atomic_sparse_kmers_count
                                            .fetch_add(1, atomic::Ordering::Relaxed);
                                    }
                                }
                                None => {
                                    // repeated part
                                    partition_kmers // separate the kmers per partition
                                        .entry(partition_index)
                                        .or_default()
                                        .push((kmer_hash, log_abundance, path_num, chunk_index));
                                    atomic_sparse_kmers_count
                                        .fetch_add(1, atomic::Ordering::Relaxed);
                                }
                            }

                            if partition_kmers.len() >= max_map_size {
                                flush_map(
                                    &mut partition_kmers,
                                    partitioned_bf_size,
                                    color_number,
                                    abundance_number,
                                    bloom_filters,
                                );
                            }
                        }
                    }
                }
                Err(e) => eprintln!("Error processing fasta: {}", e),
            }
        }
        // Flush remaining k-mers in the map by calling the earlier closure
        flush_map(
            &mut partition_kmers,
            partitioned_bf_size,
            color_number,
            abundance_number,
            bloom_filters,
        );
    })?;
    // flush the dense indexes from sparse k-mers after each file *in the first chunk*
    if chunk_index == 0 {
        if let Some(dense_indexes) = maybe_dense_indexes {
            let mut partition_kmers: HashMap<usize, Vec<(u64, u16, usize, usize)>> = HashMap::new(); // keep kmer info to fill BFs
            for (partition_index, hashmap) in dense_indexes.iter().enumerate() {
                let mut kmers_to_remove: Vec<u64> = Vec::new();
                let mut dense_index = hashmap // select the correct BF for the given partition
                    .lock()
                    .expect("Failed to lock bloom filter");
                for (kmer_hash, abundance_vector) in dense_index.iter() {
                    let number_of_zeros = count_zeros(abundance_vector, path_num)
                        .expect("Unexpected behaviour in the dense index access");
                    if number_of_zeros > threshold {
                        let color_count = path_num - number_of_zeros + 1;
                        atomic_dense_kmers_count
                            .fetch_sub(color_count as u64, atomic::Ordering::Relaxed);
                        atomic_sparse_kmers_count
                            .fetch_add(color_count as u64, atomic::Ordering::Relaxed);
                        kmers_to_remove.push(*kmer_hash);
                        for (path_index, log_abundance) in
                            abundance_vector.iter().take(path_num).enumerate()
                        {
                            if *log_abundance > 0 {
                                partition_kmers // separate the kmers per partition
                                    .entry(partition_index)
                                    .or_default()
                                    .push((
                                        *kmer_hash,
                                        (*log_abundance - 1) as u16,
                                        path_index,
                                        chunk_index,
                                    ));
                            }
                        }
                        if partition_kmers.len() >= max_map_size {
                            flush_map(
                                &mut partition_kmers,
                                partitioned_bf_size,
                                color_number,
                                abundance_number,
                                bloom_filters,
                            );
                        }
                    }
                }
                kmers_to_remove.into_iter().for_each(|key| {
                    dense_index.remove(&key);
                });
                dense_index.shrink_to_fit();
            }
            // Flush remaining k-mers in the map by calling the earlier closure
            flush_map(
                &mut partition_kmers,
                partitioned_bf_size,
                color_number,
                abundance_number,
                bloom_filters,
            );
        }
    }
    kmer_counts_vector
        .lock()
        .expect("Failed to lock the counts vector")[path_num_global] = kmer_count; // add the kmer count
    Ok(())
}

pub fn process_fasta_record(
    result: Result<fasta::Record, io::Error>,
    base: f64,
    abundance_max: u16,
    header_type: &HeaderType,
) -> Result<(Vec<u8>, u16, u16), io::Error> {
    let record = result.expect("error during fasta parsing");
    let header_option = record.desc();
    let header = header_option.unwrap_or("no header found");
    let count_value = match extract_count(header, header_type) {
        Ok(count_value) => count_value,
        Err(_) => return Err(io::Error::other("count not found!")),
    };

    // compute the lossy abundance value
    let log_abundance = if count_value > 0 {
        compute_log_abundance(count_value, base, abundance_max)
    } else {
        666
    };
    let seq = record.seq().to_vec();
    Ok((seq, log_abundance, count_value))
}

fn count_zeros(abundance_vector: &[u8], max_index: usize) -> io::Result<usize> {
    Ok(abundance_vector
        .iter()
        .take(max_index + 1)
        .filter(|&&val| val == 0)
        .count())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_process_fasta_record() {
        let fasta_input =
            ">0 ka:f:3.3   L:+:5:+ L:+:5392806:+  L:-:1:-\nAGGAGTAGATACCAGAGATAACGATACAGGTGCGA\n";

        let reader = fasta::Reader::new(std::io::Cursor::new(fasta_input));
        let result = reader.records().next().expect("failed to read record");

        let base = 2.0;
        let processed = process_fasta_record(result, base, 65535, &HeaderType::Logan);

        assert!(processed.is_ok(), "processing failed");
        let (seq, log_abundance, _) = processed.unwrap();
        let expected_seq = b"AGGAGTAGATACCAGAGATAACGATACAGGTGCGA".to_vec();
        let expected_log_abundance = 1; // log2(3) rounds to 1

        assert_eq!(seq, expected_seq, "sequence mismatch");
        assert_eq!(
            log_abundance, expected_log_abundance,
            "log abundance mismatch"
        );
    }
}
