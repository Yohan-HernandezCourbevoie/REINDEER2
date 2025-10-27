use bio::io::fasta;
use std::collections::HashMap;
use std::io;
use std::sync::atomic::Ordering;
use std::sync::{atomic, Arc, Mutex};

use crate::reindeer2::{
    compute_log_abundance, extract_count, kmer_minimizers_seq_level, process_fasta_in_batches,
    read_file, HeaderType,
};
use crate::reindeer2::{dense_index::DenseIndex, filter::Filters};

// --- INDEX FUNCTIONS ---

pub fn process_fasta_file(
    path: &str,
    maybe_dense_indexes: &Option<Arc<DenseIndex>>,
    bloom_filters: &Arc<Filters>,
    k: usize,
    m: usize,
    partition_number: usize,
    color_number_global: usize,
    threshold: usize,
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

    process_fasta_in_batches(reader, 10_000, |batch| {
        // read file 10_000 lines at once
        let mut partition_kmers: HashMap<usize, Vec<(u64, u16, usize, usize)>> = HashMap::new(); // keep kmer info to fill BFs

        // this part reads the batches of sequences and records kmers info until the structure is too large in memory
        for record in batch {
            let processed = process_fasta_record(record, base, abundance_max, &header_type); //read fasta
            match processed {
                Ok((seq, log_abundance, count_value)) => {
                    if log_abundance != 666 {
                        // case where the abundance value of the kmers in the unitigs file was < 1
                        atomic_record_count.fetch_add(1, Ordering::Relaxed);
                        let seq_str = std::str::from_utf8(&seq).expect("Invalid UTF-8 sequence");

                        for (kmer_hash, minimizer) in
                            kmer_minimizers_seq_level(seq_str.as_bytes(), k, m, canonical)
                        {
                            kmer_count += count_value as usize;
                            //for (kmer_hash, (minimizer, _)) in nt_hash_iterator.zip(min_iter) { // iterate on both minimizer and hash for each kmer
                            let partition_index = (minimizer % (partition_number as u64)) as usize;
                            total_kmers.fetch_add(1, Ordering::Relaxed);

                            // write in the dense index if the k-mer can be dense, else, put it in the hashmap for sparses
                            let inserted = match maybe_dense_indexes {
                                Some(dense_indexes) => dense_indexes.insert_if_dense(
                                    partition_index,
                                    kmer_hash,
                                    path_num_global,
                                    threshold,
                                    log_abundance,
                                    color_number_global,
                                ),
                                None => false,
                            };
                            if inserted {
                                atomic_dense_kmers_count.fetch_add(1, Ordering::Relaxed);
                            } else {
                                // write the k-mer in a file of sparse k-mer from this color
                                partition_kmers.entry(partition_index).or_default().push((
                                    kmer_hash,
                                    log_abundance,
                                    path_num,
                                    chunk_index,
                                ));
                                atomic_sparse_kmers_count.fetch_add(1, Ordering::Relaxed);
                            }

                            if partition_kmers.len() >= max_map_size {
                                bloom_filters
                                    .extend_by_draining_partitions_map(&mut partition_kmers);
                            }
                        }
                    }
                }
                Err(e) => eprintln!("Error processing fasta: {}", e),
            }
        }
        // Flush remaining k-mers in the map by calling the earlier closure
        bloom_filters.extend_by_draining_partitions_map(&mut partition_kmers);
    })?;
    // flush the dense indexes from sparse k-mers after each file *in the first chunk*
    if chunk_index == 0 {
        if let Some(dense_indexes) = maybe_dense_indexes {
            dense_indexes.truc(
                bloom_filters,
                path_num,
                threshold,
                max_map_size,
                atomic_dense_kmers_count,
                atomic_sparse_kmers_count,
                chunk_index,
            );
        }
    }
    kmer_counts_vector
        .lock()
        .expect("Failed to lock the counts vector")[path_num_global] = kmer_count; // add the kmer count
    Ok(())
}

pub fn process_fasta_record(
    record: &fasta::Record,
    base: f64,
    abundance_max: u16,
    header_type: &HeaderType,
) -> Result<(Vec<u8>, u16, u16), io::Error> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_process_fasta_record() {
        let fasta_input =
            ">0 ka:f:3.3   L:+:5:+ L:+:5392806:+  L:-:1:-\nAGGAGTAGATACCAGAGATAACGATACAGGTGCGA\n";

        let reader = fasta::Reader::new(std::io::Cursor::new(fasta_input));
        let result = reader.records().next().expect("failed to read record");
        let result = result.expect("error during fasta parsing");

        let base = 2.0;
        let processed = process_fasta_record(&result, base, 65535, &HeaderType::Logan);

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
