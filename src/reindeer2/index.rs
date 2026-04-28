use bio::io::fasta;
use log::warn;
use std::collections::HashMap;
use std::io;
use std::num::NonZero;
use std::sync::atomic::Ordering;
use std::sync::{Arc, Mutex, atomic};

use crate::reindeer2::{
    HeaderType, compute_log_abundance, extract_count,
    minimizer_iter::{KmerMinimizerIteratorError, Sampler, kmer_minimizers_sampled},
    process_fasta_in_batches, read_file,
};
use crate::reindeer2::{dense_index::DenseIndex, storage::filters::Filters};

// --- INDEX FUNCTIONS ---

pub fn process_fasta_file<S>(
    path: &str,
    maybe_dense_indexes: &Option<Arc<DenseIndex>>,
    bloom_filters: &Filters,
    k: usize,
    z: usize,
    m: usize,
    partition_number: usize,
    color_number_global: usize,
    threshold: usize,
    abundance_min: u16,
    abundance_max: NonZero<u16>,
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
    count_right_after_angle_bracket: bool,
    sampler: &S,
) -> io::Result<()>
where
    S: Sampler,
{
    let atomic_record_count = atomic::AtomicU64::new(0);
    let mut kmer_count: usize = 0;
    let reader = read_file(path)?;

    process_fasta_in_batches(reader, 10_000, |batch| {
        // read file 10_000 lines at once
        let mut partition_kmers: HashMap<usize, Vec<(u64, u16, usize, usize)>> = HashMap::new(); // keep kmer info to fill BFs

        // this part reads the batches of sequences and records kmers info until the structure is too large in memory
        for record in batch {
            let processed = process_fasta_record(
                record,
                base,
                abundance_min,
                abundance_max,
                &header_type,
                count_right_after_angle_bracket,
            ); //read fasta
            match processed {
                Ok((seq, log_abundance, count_value)) => {
                    // TODO discuss we should at least print a warning here if it is None
                    if let Some(log_abundance) = log_abundance {
                        // case where the abundance value of the kmers in the unitigs file was < 1
                        atomic_record_count.fetch_add(1, Ordering::Relaxed);
                        let seq_str = std::str::from_utf8(&seq).expect("Invalid UTF-8 sequence");

                        let smer_minimizers = match kmer_minimizers_sampled(
                            seq_str.as_bytes(),
                            k - z, // use k-z mers instead of k-mers
                            m,
                            canonical,
                            sampler,
                        ) {
                            Ok(iterator) => iterator,
                            Err(KmerMinimizerIteratorError::SequenceTooSmall { k }) => {
                                eprintln!(
                                    "Warning: when indexing file {path}, the read {seq_str} will be ignored. To be indexed, its length (={}) must be greater or equal to k (={})",
                                    seq_str.len(),
                                    k
                                );
                                continue;
                            }
                        };

                        kmer_count += count_value as usize * (seq.len() - k);
                        for (smer_hash, minimizer) in smer_minimizers {
                            let partition_index = (minimizer % (partition_number as u64)) as usize;
                            total_kmers.fetch_add(1, Ordering::Relaxed);

                            // write in the dense index if the k-mer can be dense, else, put it in the hashmap for sparses
                            let inserted = match maybe_dense_indexes {
                                Some(dense_indexes) => dense_indexes.insert_if_dense(
                                    partition_index,
                                    smer_hash,
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
                                    smer_hash,
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
                Err(e) => {
                    warn!("Error processing fasta: {}", e);
                    eprintln!("Error processing fasta: {}", e)
                }
            }
        }
        // Flush remaining k-mers in the map by calling the earlier closure
        bloom_filters.extend_by_draining_partitions_map(&mut partition_kmers);
    })?;
    // flush the dense indexes from sparse k-mers after each file *in the first chunk*
    if chunk_index == 0
        && let Some(dense_indexes) = maybe_dense_indexes
    {
        dense_indexes.remove_sparse_entries(
            bloom_filters,
            path_num,
            threshold,
            max_map_size,
            atomic_dense_kmers_count,
            atomic_sparse_kmers_count,
            chunk_index,
        );
    }
    kmer_counts_vector
        .lock()
        .expect("Failed to lock the counts vector")[path_num_global] = kmer_count; // add the kmer count
    Ok(())
}

pub fn process_fasta_record(
    record: &fasta::Record,
    base: f64,
    abundance_min: u16,
    abundance_max: NonZero<u16>,
    header_type: &HeaderType,
    count_right_after: bool,
) -> Result<(Vec<u8>, Option<u16>, u16), io::Error> {
    let header = if count_right_after {
        record.id()
    } else {
        match record.desc() {
            Some(header) => header,
            None => return Err(io::Error::other("no header found")),
        }
    };
    let count_value = match extract_count(header, header_type) {
        Ok(count_value) => count_value,
        Err(_) => return Err(io::Error::other("count not found!")),
    };

    // compute the lossy abundance value
    let abundance = if count_value > abundance_min {
        NonZero::new(count_value)
    } else {
        None
    };
    let log_abundance = abundance.map(|value| compute_log_abundance(value, base, abundance_max));
    let seq = record.seq().to_vec();
    Ok((seq, log_abundance, count_value))
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use super::*;

    // TODO test abundance min/max
    #[test]
    fn test_process_fasta_record() {
        let fasta_input =
            ">0 ka:f:3.3   L:+:5:+ L:+:5392806:+  L:-:1:-\nAGGAGTAGATACCAGAGATAACGATACAGGTGCGA\n";

        let reader = fasta::Reader::new(std::io::Cursor::new(fasta_input));
        let result = reader.records().next().expect("failed to read record");
        let result = result.expect("error during fasta parsing");

        let base = 2.0;
        let max = NonZero::new(65535).unwrap();
        // TODO add a test using true
        let processed = process_fasta_record(&result, base, 0, max, &HeaderType::Logan, false);

        assert!(processed.is_ok(), "processing failed");
        let (seq, log_abundance, _) = processed.unwrap();
        let expected_seq = b"AGGAGTAGATACCAGAGATAACGATACAGGTGCGA".to_vec();
        let expected_log_abundance = 1; // log2(3) rounds to 1

        assert_eq!(seq, expected_seq, "sequence mismatch");
        assert_eq!(
            log_abundance.unwrap(),
            expected_log_abundance,
            "log abundance mismatch"
        );
    }
}
