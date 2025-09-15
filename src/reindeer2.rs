use std::collections::{HashMap, VecDeque};
use std::io::{self, BufRead, BufReader,Write, BufWriter,Read};
use std::panic;
use std::fs::{self, File};
use std::path::Path;
use std::time::Instant;
use std::sync::{Arc, Mutex, atomic};
use std::error::Error;

use flate2::read::GzDecoder;
use rayon::prelude::*;
use bio::io::fasta;
use zstd::stream::decode_all;
use roaring::RoaringBitmap;
use nthash::NtHashIterator;
use csv::Writer;
use thousands::Separable;
use num_format::{Locale, ToFormattedString};



// === INDEX ===

// --- MAIN FUNCTION ---

pub fn build_index(
    file_paths: Vec<String>, 
    kmer_size: usize, 
    minimizer_size: usize, 
    bf_size: u64,
    partition_number: usize,
    color_nb: usize,
    abundance_number: usize,
    abundance_max: u16,
    output_dir: &str,
    dense_option: bool,
    threshold: usize,
    debug: bool,
) -> io::Result<(Vec<String>, String)> {
    let mut total_kmers = atomic::AtomicU64::new(0);
    let mut atomic_dense_kmers_count = atomic::AtomicU64::new(0);
    let mut atomic_sparse_kmers_count = atomic::AtomicU64::new(0);
    let mut atomic_sparse_one_seen = atomic::AtomicU64::new(0);
    let mut atomic_sparse_fp_seen = atomic::AtomicU64::new(0);
    let kmer_counts_vector: Arc<Mutex<Vec<usize>>> = Arc::new(Mutex::new(vec![0; color_nb]));
    let (chunks, color_chunks) = split_fof(&file_paths)?;
    let base = compute_base(abundance_number, abundance_max);
    if debug {
        println!("Using log base {}", base);
    }
    let partitioned_bf_size = (bf_size as usize) / partition_number;
    if debug {
        println!("In debug mode... the tool may take (much) longer than usual.");
    }
    println!("Initializing Bloom filter slices...");

    let (_, dir_path) = create_dir_and_files(partition_number, output_dir)?;

    // Shared data structures protected by Mutex for safe parallel access
    let maybe_dense_indexes: Option<Arc<Vec<Mutex<HashMap<u64, Vec<u8>>>>>> = 
        if dense_option {
            Some(Arc::new(
            (0..partition_number)
                .map(|_| Mutex::new(HashMap::with_capacity(200_000_000/partition_number)))
                .collect::<Vec<_>>()
            ))
        } else {
            None
        };

    let bloom_filters = Arc::new(
        (0..partition_number)
            .map(|_| Mutex::new(RoaringBitmap::new()))
            .collect::<Vec<_>>()
    );

    for (chunk_i, chunk) in chunks.iter().enumerate() {
        // For each file in this chunk, process in *parallel* (soon)
        // TODO build the appropriate iterator to parallelize (or not) if dense is set
        chunk.par_iter().enumerate().for_each(|(path_num, path)| {
            match std::fs::metadata(path) {
                Ok(metadata) => {
                    if metadata.is_file() {
                        let reader = match read_file(path) {
                            Ok(r) => r,
                            Err(e) => {
                                eprintln!("Failed to open file {}: {}", path, e);
                                return;
                            }
                        };
                        let fasta_reader = fasta::Reader::new(reader);
                        let first_record = fasta_reader.records().next();

                        if let Some(Ok(record)) = first_record {
                            let header_option = record.desc();
                            let header = header_option.unwrap_or("no header found");
                            let h_type = match determine_header_type(header) {
                                Ok(ht) => ht,
                                Err(e) => {
                                    eprintln!("Unsupported header type ({}): {}", path, e);
                                    return;
                                }
                            };

                            if let Err(e) = process_fasta_file(
                                path,
                                &maybe_dense_indexes,
                                &bloom_filters,
                                kmer_size,
                                minimizer_size,
                                partitioned_bf_size,
                                partition_number,
                                color_chunks[chunk_i],
                                color_nb,
                                threshold,
                                abundance_number,
                                abundance_max,
                                path_num % color_chunks[0],
                                path_num,
                                base,
                                chunk_i,
                                h_type, 
                                1_000_000, // max size for flushing k-mers to bloom filter
                                &total_kmers,
                                &atomic_dense_kmers_count,
                                &atomic_sparse_kmers_count,
                                &kmer_counts_vector,
                            ) {
                                eprintln!("Error processing {}: {}", path, e);
                            }
                        } else {
                            eprintln!("Failed to determine header type for {}", path);
                        }
                    } else {
                        eprintln!("Path {} exists but is not a file", path);
                    }
                }
                Err(_) => eprintln!("Path {} does not exist", path),
            }
        });
        if debug {
            update_sparse_counts(&bloom_filters, &atomic_sparse_one_seen, &atomic_sparse_fp_seen, abundance_number)?;
        }
        if let Err(e) = write_bloom_filters_to_disk(
            &dir_path, 
            &bloom_filters,
            &color_chunks,
            partition_number,
            chunk_i
        ) {
            eprintln!("Error writing Bloom filters for chunk {}: {}", chunk_i, e);
        }
        println!("Chunk {} done", chunk_i);
    };

    // After processing all chunks, write the dense indexes to disk
    if let Some(dense_indexes) = maybe_dense_indexes {
        write_dense_indexes_to_disk(&dir_path, &dense_indexes)?;
    }

    write_kmer_counts_to_disk(&dir_path, &kmer_counts_vector)?;

    if debug {
        // k-mers repartition between dense and sparse index
        println!("The index contains {:?} 'dense' k-mers and {:?} 'sparse' k-mers (total k-mers: {:?})", 
            atomic_dense_kmers_count.get_mut().separate_with_commas(), 
            atomic_sparse_kmers_count.get_mut().separate_with_commas(), 
            total_kmers.get_mut().separate_with_commas());

        let ones = atomic_sparse_one_seen.get_mut();
        let silent = *atomic_sparse_kmers_count.get_mut() - *ones;
        let fp = atomic_sparse_fp_seen.get_mut();
        // k-mers indexed in the sparse index, FP silent and FP seen
        println!("Among the {:?} k-mers added in the 'sparse' index, {:?} encountered hash collisions ({:?} silent and {:?} misleading).", 
            atomic_sparse_kmers_count.get_mut().separate_with_commas(), 
            (silent+*fp).separate_with_commas(),
            silent.separate_with_commas(), 
            fp.separate_with_commas());
    }

    // Merge the chunk-based bloom filters, if needed
    if chunks.len() > 1 {
        merge_all_partitions(
            &dir_path,
            &dir_path, // output
            partitioned_bf_size,
            abundance_number,
            color_chunks,
            partition_number,
            color_nb,
        )?;
    } else {
        // If there was only one chunk, rename files directly
        for partition_idx in 0..partition_number {
            let input_path = format!("{}/partition_bloom_filters_c0_p{}.txt", dir_path, partition_idx);
            let output_path = format!("{}/partition_bloom_filters_p{}.txt", dir_path, partition_idx);
            std::fs::rename(&input_path, &output_path)?;
        }
    } 

    // write partition info to a CSV or your desired format
    let _ = write_partition_to_csv(
        &dir_path,
        kmer_size,
        minimizer_size,
        bf_size,
        partition_number,
        color_nb,
        abundance_number,
        abundance_max,
        dense_option,
    );

    Ok((file_paths, dir_path))
}

pub fn build_index_muset(
    unitigs_file: String, 
    matrix_file: String, 
    kmer_size: usize, 
    minimizer_size: usize, 
    bf_size: u64,
    partition_number: usize,
    color_nb: usize,
    abundance_number: usize,
    abundance_max: u16,
    output_dir: &str,
    dense_option: bool,
    threshold: usize,
    debug: bool,
) -> io::Result<(Vec<String>, String)> {
    let mut total_kmers = atomic::AtomicU64::new(0);
    let mut atomic_dense_kmers_count = atomic::AtomicU64::new(0);
    let mut atomic_sparse_kmers_count = atomic::AtomicU64::new(0);
    let mut atomic_sparse_one_seen = atomic::AtomicU64::new(0);
    let mut atomic_sparse_fp_seen = atomic::AtomicU64::new(0);
    let kmer_counts_vector: Arc<Mutex<Vec<usize>>> = Arc::new(Mutex::new(vec![color_nb, 0]));
    let base = compute_base(abundance_number, abundance_max);
    let max_map_size = 1_000_000;
    if debug {
        println!("Using log base {}", base);
    }
    let partitioned_bf_size = (bf_size as usize) / partition_number;
    if debug {
        println!("In debug mode... the tool may take (much) longer than usual.");
    }
    println!("Initializing Bloom filter slices...");

    let (_, dir_path) = create_dir_and_files(partition_number, output_dir)?;

    // Shared data structures protected by Mutex for safe parallel access
    let maybe_dense_indexes: Option<Arc<Vec<Mutex<HashMap<u64, Vec<u8>>>>>> = 
        if dense_option {
            Some(Arc::new(
            (0..partition_number)
                .map(|_| Mutex::new(HashMap::with_capacity(200_000_000/partition_number)))
                .collect::<Vec<_>>()
            ))
        } else {
            None
        };

    let bloom_filters = Arc::new(
        (0..partition_number)
            .map(|_| Mutex::new(RoaringBitmap::new()))
            .collect::<Vec<_>>()
    );

    let unitigs_reader = match read_file(&unitigs_file) {
        Ok(r) => r,
        Err(e) => {
            panic!("Failed to open file {}: {}", unitigs_file, e);
        }
    };
    let matrix_reader = match read_file(&matrix_file) {
        Ok(r) => r,
        Err(e) => {
            panic!("Failed to open file {}: {}", matrix_file, e);
        }
    };
    let len_matrix = matrix_reader.lines().count();
    let len_unitigs = unitigs_reader.lines().count();
    if len_matrix * 2 != len_unitigs {
        panic!("The number of unitigs doesn't match the matrix size.");
    }

    let mut matrix_lines = read_file(&matrix_file)?.lines();
    let mut unitigs_lines = read_file(&unitigs_file)?.lines();

    let mut partition_kmers: HashMap<usize, Vec<(u64, u16, usize, usize)>> = HashMap::new(); // keep kmer info to fill BFs

    let mut kmer_counts_vector_locked = kmer_counts_vector.lock().expect("Failed to lock the kmer counter hashmap");
    for _ in 0..len_matrix {
        let abundance_line = matrix_lines.next().unwrap().unwrap();
        let mut abundance_iter = abundance_line.trim().split(" ").into_iter();
        let unitig_id: &str = abundance_iter.next().unwrap();

        let abundance_vector_plusone = abundance_iter.enumerate().map(|(i, ab_value_str)| {
            let ab_value = ab_value_str.parse::<f64>().expect("Failed to parse the abundance values in the matrix") as u16;
            if ab_value == 0 {
                0
            } else {
                kmer_counts_vector_locked[i] += 1; // select the count corresponding to the right color
                (compute_log_abundance(
                    ab_value,
                    base,
                    abundance_max) + 1) as u8
            }})
            .collect::<Vec<u8>>();

        let tmp_abundance_vector = abundance_vector_plusone.clone();
        let number_of_zeros: usize = tmp_abundance_vector.into_iter().filter(|x| *x == 0).count();

        let unitig_line = unitigs_lines.next().unwrap().unwrap();
        if &unitig_line.split(" ").next().unwrap()[1..] != unitig_id {
            println!("{} vs {}", &unitig_line.split(" ").next().unwrap()[1..], unitig_id);
            panic!("The unitig id dooesn't match between the unitigs file and the matrix file.")
        }

        let unitig_line = unitigs_lines.next().unwrap().unwrap();
        let seq_str = unitig_line.trim();
        for (kmer_hash, minimizer) in kmer_minimizers_seq_level(seq_str.as_bytes(), kmer_size, minimizer_size) { 
            let partition_index = (minimizer % (partition_number as u64)) as usize;
            let new_kmers_added = (abundance_vector_plusone.len() - number_of_zeros) as u64;
            total_kmers.fetch_add(new_kmers_added, atomic::Ordering::Relaxed);
            
            if dense_option && number_of_zeros <= threshold {
                match &maybe_dense_indexes {
                    Some(dense_indexes) => {
                        let mut dense_index = dense_indexes[partition_index] 
                            .lock()
                            .expect("Failed to lock the dense index");
                        atomic_dense_kmers_count.fetch_add(new_kmers_added, atomic::Ordering::Relaxed);
                        dense_index.insert(kmer_hash, abundance_vector_plusone.clone());
                    },
                    None => panic!("Failed to build the dense index.")
                }
            } else {
                atomic_sparse_kmers_count.fetch_add(new_kmers_added, atomic::Ordering::Relaxed);
                let kmers_entry = partition_kmers // separate the kmers per partition
                    .entry(partition_index)
                    .or_default();

                for (path_num, log_plusone) in abundance_vector_plusone.iter().enumerate() {
                    let real_log_abundance = (log_plusone - 1) as u16;
                    if real_log_abundance > 0 {
                        kmers_entry.push((
                            kmer_hash,
                            real_log_abundance,
                            path_num,
                            0 // chunk index
                        ))
                    }
                }

            }
            if partition_kmers.len() >= max_map_size {
                flush_map_into_bfs(&mut partition_kmers, &bloom_filters, partitioned_bf_size, color_nb, abundance_number)?;
            }
        }
    }
    // Index remaining kmers
    flush_map_into_bfs(&mut partition_kmers, &bloom_filters, partitioned_bf_size, color_nb, abundance_number)?;
        

        


        

    if debug {
        update_sparse_counts(&bloom_filters, &atomic_sparse_one_seen, &atomic_sparse_fp_seen, abundance_number)?;
    }
    if let Err(e) = write_bloom_filters_to_disk(
        &dir_path, 
        &bloom_filters,
        &vec![color_nb],
        partition_number,
        0
    ) {
        eprintln!("Error writing Bloom filters on disk: {}", e);
    }

    // After processing all chunks, write the dense indexes to disk
    if let Some(dense_indexes) = maybe_dense_indexes {
        write_dense_indexes_to_disk(&dir_path, &dense_indexes)?;
    }

    write_kmer_counts_to_disk(&dir_path, &kmer_counts_vector)?;

    if debug {
        // k-mers repartition between dense and sparse index
        println!("The index contains {:?} 'dense' k-mers and {:?} 'sparse' k-mers (total k-mers: {:?})", 
            atomic_dense_kmers_count.get_mut().separate_with_commas(), 
            atomic_sparse_kmers_count.get_mut().separate_with_commas(), 
            total_kmers.get_mut().separate_with_commas());

        let ones = atomic_sparse_one_seen.get_mut();
        let silent = *atomic_sparse_kmers_count.get_mut() - *ones;
        let fp = atomic_sparse_fp_seen.get_mut();
        // k-mers indexed in the sparse index, FP silent and FP seen
        println!("Among the {:?} k-mers added in the 'sparse' index, {:?} encountered hash collisions ({:?} silent and {:?} misleading).", 
            atomic_sparse_kmers_count.get_mut().separate_with_commas(), 
            (silent+*fp).separate_with_commas(),
            silent.separate_with_commas(), 
            fp.separate_with_commas());
    }
    

    // Merge the chunk-based bloom filters, if needed
    // If there was only one chunk, rename files directly
    for partition_idx in 0..partition_number {
        let input_path = format!("{}/partition_bloom_filters_c0_p{}.txt", dir_path, partition_idx);
        let output_path = format!("{}/partition_bloom_filters_p{}.txt", dir_path, partition_idx);
        std::fs::rename(&input_path, &output_path)?;
    }

    // write partition info to a CSV or your desired format
    let _ = write_partition_to_csv(
        &dir_path,
        kmer_size,
        minimizer_size,
        bf_size,
        partition_number,
        color_nb,
        abundance_number,
        abundance_max,
        dense_option,
    );

    Ok((vec!["".to_string()], dir_path))
}


// --- INDEX FUNCTIONS ---

fn process_fasta_file(
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
) -> io::Result<()> {
    let atomic_record_count = atomic::AtomicU64::new(0);
    let mut kmer_count: usize = 0;
    let reader = read_file(path)?;

    // this part fills the BFs per partition
    fn flush_map(partition_kmers: &mut HashMap<usize, Vec<(u64, u16, usize, usize)>>, partitioned_bf_size: usize, color_number: usize, abundance_number: usize, bloom_filters: &Arc<Vec<Mutex<RoaringBitmap>>>) {
        for (partition_index, kmers) in partition_kmers.drain() {// iterates and empties the hash map when needed
            

            let mut kmer_hashes_to_update = Vec::with_capacity(kmers.len());
            for (kmer_hash, log_abundance, path_num, _chunk_index) in kmers { // select the bit in the BF
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

    process_fasta_in_batches(reader, 10_000, |batch| { // read file 10_000 lines at once
        let mut partition_kmers: HashMap<usize, Vec<(u64, u16, usize, usize)>> = HashMap::new(); // keep kmer info to fill BFs

        // this part reads the batches of sequences and records kmers info until the structure is too large in memory
        for record in batch {
            let processed = process_fasta_record(Ok(record), base, abundance_max, &header_type); //read fasta
            match processed {
                Ok((seq, log_abundance, count_value)) => {
                    if log_abundance != 666 { // case where the abundance value of the kmers in the unitigs file was < 1
                        atomic_record_count.fetch_add(1, atomic::Ordering::Relaxed);
                        let seq_str = std::str::from_utf8(&seq).expect("Invalid UTF-8 sequence");

                        for (kmer_hash, minimizer) in kmer_minimizers_seq_level(seq_str.as_bytes(), k, m) {
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
                                        if let Some(abundance_vector) = dense_index.get_mut(&kmer_hash) {
                                            abundance_vector[path_num_global] = (log_abundance + 1) as u8;
                                        }
                                        atomic_dense_kmers_count.fetch_add(1, atomic::Ordering::Relaxed);
                                    } else if path_num_global <= threshold {
                                        // create a new abundance vector for the k-mer
                                        let mut abundance_vector: Vec<u8> = Vec::with_capacity(color_number_global);
                                        abundance_vector.resize(color_number_global, 0);
                                        abundance_vector[path_num_global] = (log_abundance + 1) as u8;
                                        dense_index.insert(kmer_hash, abundance_vector);
                                        atomic_dense_kmers_count.fetch_add(1, atomic::Ordering::Relaxed);
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
                                        atomic_sparse_kmers_count.fetch_add(1, atomic::Ordering::Relaxed);
                                    }
                                },
                                None => {
                                    // repeated part
                                    partition_kmers // separate the kmers per partition
                                        .entry(partition_index)
                                        .or_default()
                                        .push((
                                            kmer_hash,
                                            log_abundance,
                                            path_num,
                                            chunk_index,
                                        ));
                                    atomic_sparse_kmers_count.fetch_add(1, atomic::Ordering::Relaxed);
                                }
                            }

                            
                            if partition_kmers.len() >= max_map_size {
                                flush_map(&mut partition_kmers, partitioned_bf_size, color_number, abundance_number, bloom_filters);
                            }
                        }
                    }
                }
                Err(e) => eprintln!("Error processing fasta: {}", e),
            }
        }
        // Flush remaining k-mers in the map by calling the earlier closure
        flush_map(&mut partition_kmers, partitioned_bf_size, color_number, abundance_number, bloom_filters);
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
                    let number_of_zeros = count_zeros(abundance_vector, path_num).expect("Unexpected behaviour in the dense index access");
                    if number_of_zeros > threshold {
                        let color_count = path_num - number_of_zeros + 1;
                        atomic_dense_kmers_count.fetch_sub(color_count as u64, atomic::Ordering::Relaxed);
                        atomic_sparse_kmers_count.fetch_add(color_count as u64, atomic::Ordering::Relaxed);
                        kmers_to_remove.push(*kmer_hash);
                        for (path_index, log_abundance) in abundance_vector.iter().take(path_num).enumerate() {
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
                            flush_map(&mut partition_kmers, partitioned_bf_size, color_number, abundance_number, bloom_filters);
                        }
                    }
                }
                kmers_to_remove.into_iter().for_each(|key| {dense_index.remove(&key);});
                dense_index.shrink_to_fit();
            }
            // Flush remaining k-mers in the map by calling the earlier closure
flush_map(&mut partition_kmers, partitioned_bf_size, color_number, abundance_number, bloom_filters);        }
    }
    kmer_counts_vector 
        .lock()
        .expect("Failed to lock the counts vector")
        [path_num_global] = kmer_count; // add the kmer count
    Ok(())
}

fn process_fasta_in_batches<R: io::BufRead>(
    reader: R,
    batch_size: usize,
    mut process_batch: impl FnMut(Vec<fasta::Record>),
) -> io::Result<()> {
    let fasta_reader = fasta::Reader::new(reader);
    let mut batch = Vec::with_capacity(batch_size);

    for result in fasta_reader.records() {
        let record = result?; 
        batch.push(record);

        if batch.len() >= batch_size {
            process_batch(batch);
            batch = Vec::with_capacity(batch_size); // reset
        }
    }

    // final batch
    if !batch.is_empty() {
        process_batch(batch);
    }

    Ok(())
}

fn process_fasta_record (
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

fn count_zeros(
    abundance_vector: &[u8],
    max_index: usize,
) -> io::Result<usize> {
    Ok(abundance_vector
        .iter()
        .take(max_index+1)
        .filter(|&&val| val == 0)
        .count())
}



// === QUERY ===

// --- MAIN FUNCTION ---

pub fn query_index(
    fasta_file: &str,
    bf_dir: &str,
    output_file: &str,
    color_graph: bool,
    normalize_option: bool,
    coverage: f32,
) -> io::Result<()> {
    //load index metadata from CSV
    let (k, m, bf_size, partition_number, color_number, abundance_number, abundance_max, dense_option) =
        read_partition_from_csv(bf_dir, "index_info.csv")?;
    println!("Loaded index metadata for query.");
   
    
    let query_results = query_sequences_in_batches(
        fasta_file,
        bf_dir,
        k,
        m,
        bf_size,
        partition_number,
        color_number,
        abundance_number,
        abundance_max,
        10_000_000,
        &output_file,
        color_graph,
        dense_option,
        normalize_option,
        coverage,
    )?;
    println!("Writing results in {}", output_file);
    //write_query_results_to_csv(&query_results, bf_dir)
    Ok(query_results)
}


// --- QUERY FUNCTIONS ---

fn query_sequences_in_batches(
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
    color_graph: bool,
    dense_option: bool,
    normalize_option: bool,
    coverage: f32,
) -> io::Result<()> {
    let reader = read_file(fasta_file)?;
    let mut writer = BufWriter::new(File::create(output_file)?);
    let base = compute_base(abundance_number, abundance_max);

    // write the headr of the output csv file
    writeln!(writer, "header,file,abundance")?;

    // Process FASTA in chunks of `batch_size` records
    process_fasta_in_batches(reader, batch_size, |batch| {
        // Build a map from partition_index -> Vec of (sequence_id, kmer_hash).
        let mut partition_kmers: HashMap<usize, Vec<(String, u64)>> = HashMap::new();

        // Build all partition-kmers upfront
        for record in batch {
            let id = record.id();
            let desc = record.desc().unwrap_or("");
            // Build the header string only once per record
            let full_header = format!(">{} {}", id, desc).trim().to_string();
            // let full_header = format!(">{}", id).trim().to_string();

            let seq_str = std::str::from_utf8(record.seq())
                .map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidData, "Invalid UTF-8 sequence")
                })
                .unwrap(); // or handle error gracefully

            // k-mer hash + minimizer
            //let nt_hash_iterator = NtHashIterator::new(seq_str.as_bytes(), k).unwrap(); 
            //let min_iter = MinimizerBuilder::<u64>::new()
            //    .minimizer_size(m)
            //    .width((k - m + 1).try_into().unwrap())
            //    .iter(seq_str.as_bytes());

            for (kmer_hash, minimizer) in kmer_minimizers_seq_level(seq_str.as_bytes(), k, m) {
                //for (kmer_hash, (minimizer, _position)) in nt_hash_iterator.zip(min_iter) {
                let partition_index = (minimizer % (partition_number as u64)) as usize;
                partition_kmers
                    .entry(partition_index)
                    .or_default()
                    .push((full_header.clone(), kmer_hash));
            }
        }

        // --- PARALLEL PHASE: process each partition's k-mers in parallel ---
        // This will store final results for *all* sequences in this batch.
        // Key: sequence header; Value: vector of vector of abundances
        let sequence_results = partition_kmers
            .into_par_iter()
            // 1) Create a local HashMap in each thread
            .fold(
                HashMap::<String, Vec<Vec<u16>>>::new,
                |mut local_results, (partition_index, kmers)| {
                    // Load the partition's Bloom filter
                    let path_bf = format!("{}/partition_bloom_filters_p{}.txt",bf_dir, partition_index);
                    let maybe_bf = load_bloom_filter(&path_bf);

                    if let Ok((bitmap, _maybe_aux_data)) = maybe_bf {
                        let hashmap: HashMap<u64, Box<Vec<u8>>> = if dense_option {
                            let path_dense_index = format!("{}/partition_dense_index_p{}.txt", bf_dir, partition_index);
                            load_dense_index(&path_dense_index).expect(&format!("Failed to load dense index for partition {}", partition_index))
                        } else {
                            HashMap::new()
                        };
                        //  For each k-mer in this partition
                        for (sequence_id, kmer_hash) in kmers {
                            let color_abundances = 
                                if hashmap.contains_key(&kmer_hash) {
                                    let log_abundance_vector = hashmap.get(&kmer_hash).expect("failed to read the hashmap");
                                    let mut color_abundances = vec![Vec::new(); color_number];
                                    for (color, log_abundance) in log_abundance_vector.iter().enumerate() {
                                            if *log_abundance == 0 {
                                                // TODO discuss weird value
                                                color_abundances[color].push(666);
                                            } else {
                                                color_abundances[color].push((*log_abundance - 1) as usize);
                                            }
                                        };
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
                                    update_color_abundances(
                                        &bitmap,
                                        base_position,
                                        color_number,
                                        abundance_number,
                                        &mut color_abundances,
                                    );
                                    color_abundances
                                };

                            // Convert log abundances to approximate integer counts
                            let approximate_counts: Vec<Vec<u16>> = color_abundances
                                .into_iter()
                                .map(|abunds_for_color| {
                                    abunds_for_color
                                        .into_iter()
                                        .map(|log_abund| if log_abund == 666 {
                                                0
                                            } else {
                                                approximate_value(log_abund, base)
                                            })
                                        .collect()
                                })
                                .collect();

                            // Accumulate results in local_results
                            let entry = local_results
                                .entry(sequence_id.clone())
                                .or_insert_with(|| vec![Vec::new(); color_number]);
                            for (color_idx, approx_values) in approximate_counts.into_iter().enumerate() {
                                entry[color_idx].push(*approx_values.iter().min().expect("An abundance vector returned empty"));
                            }
                        }
                        
                    } else {
                        eprintln!("Failed to load Bloom filter for partition {}", partition_index);
                    }

                    local_results
                }
            )
            // 2) Reduce all local HashMaps into a single HashMap
            .reduce(
                HashMap::<String, Vec<Vec<u16>>>::new,
                merge_results,
            );

        // Now `sequence_results` has the combined data for this batch.
        // Either color a graph or compute medians and output them.
        if color_graph {
            // If your graph coloring wants to read from `sequence_results`:
            // Flush the writer to separate batch outputs if needed
            let _ = writer.flush();
            let _ = graph_coloring(
                fasta_file,
                batch_size,
                output_file,
                &sequence_results
            );
        } else {
            let kmer_counts = if normalize_option {
                load_kmer_counts_vector(bf_dir).expect("Failed to load from disk the kmer counts vector")
            } else {
                vec![color_number, 0]
            };
            // Compute medians for each sequence and each color, then write them out
            for (seq_header, color_vectors) in &sequence_results {
                for (color_idx, abund_values) in color_vectors.iter().enumerate() {
                    if !abund_values.is_empty() {
                        let mut zeros_count = 0;
                        let mut non_zero_values: Vec<u16> = Vec::new();
                        abund_values.iter().for_each(|value|
                            if *value == 0 {
                                zeros_count += 1;
                            } else {
                                non_zero_values.push(*value);
                            }
                        );
                        if !non_zero_values.is_empty() && (((zeros_count as f32) / (abund_values.len() as f32)) < coverage) {
                            let mut abund_sorted = non_zero_values.clone();
                            abund_sorted.sort_unstable();
                            let median = 
                                if abund_sorted.iter().all(|&x| x == 0) {
                                    0
                                } else if abund_sorted.len() == 1 {
                                    abund_sorted[0]
                                } else {
                                    let mid = abund_sorted.len() / 2;
                                    if abund_sorted.len() % 2 == 1 {
                                        abund_sorted[mid]
                                    } else {
                                        (abund_sorted[mid - 1] + abund_sorted[mid]) / 2
                                    }
                                };
                            if median > 0 {
                                let median = if normalize_option {
                                    median as f64 / kmer_counts[color_idx] as f64 * 1_000_000 as f64
                                } else {
                                    median as f64
                                };
                                let _ = writeln!(
                                    writer,
                                    "{},{},{}",
                                    seq_header,
                                    color_idx,
                                    median
                                );
                            }
                        }
                    }
                }
            }
            let _ = writer.flush();
        }
    })?;
    
    Ok(())
}

fn merge_results(
    mut acc: HashMap<String, Vec<Vec<u16>>>,
    local: HashMap<String, Vec<Vec<u16>>>
) -> HashMap<String, Vec<Vec<u16>>> {
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

// rewrites a bcalm-like graph so that headers have abund info (one of the possible query operations)
pub fn graph_coloring(
    fasta_file: &str,
    batch_size: usize,
    output_file: &str,
    sequence_results: &HashMap<String, Vec<Vec<u16>>>,
) -> Result<(), Box<dyn Error>> {
    let reader = read_file(fasta_file)?; // your existing read_file
    let mut writer = BufWriter::new(File::create(output_file)?);

    process_fasta_in_batches(reader, batch_size, |batch| {
        for record in batch {
            let id = record.id();
            let desc = record.desc().unwrap_or("");
            let full_header = format!(">{} {}", id, desc).trim().to_string();
            let seq_str = std::str::from_utf8(record.seq())
                .expect("Invalid UTF-8 sequence");

            // If the sequence is in sequence_results, we fetch the vec of vec
            if let Some(color_vectors) = sequence_results.get(&full_header) {
                // color_vectors is a Vec<Vec<u16>>. Each index = a color,
                // each inner Vec<u16> = all abundance values for that color
                // if no data, just write the original header
                if color_vectors.iter().all(|vals| vals.is_empty()) {
                    writeln!(writer, "{}", full_header).ok();
                    writeln!(writer, "{}", seq_str).ok();
                    continue;
                }
                // otherwise, build an augmented header
                let mut header_parts = Vec::with_capacity(color_vectors.len() + 1);
                header_parts.push(full_header.clone());

                // for each color, we do the median of all values:
                for (color_idx, vals) in color_vectors.iter().enumerate() {
                    if vals.is_empty() {
                        // skip color if it has no data
                        continue;
                    }
                    // compute median
                    let mut abund_sorted = vals.clone();
                        abund_sorted.sort_unstable();
                        let median = 
                            if abund_sorted.iter().all(|&x| x == 0) {
                                0
                            } else if abund_sorted.len() == 1 {
                                abund_sorted[0]
                            } else {
                                let mid = abund_sorted.len() / 2;
                                if abund_sorted.len() % 2 == 1 {
                                    abund_sorted[mid]
                                } else {
                                    (abund_sorted[mid - 1] + abund_sorted[mid]) / 2
                                }
                            };

                    // push e.g. "col:1:12"
                    header_parts.push(format!("col:{}:{}", color_idx, median));
                }

                // join info like
                // ">seq1 col:0:12 col:1:29"
                let new_header = header_parts.join(" ");
                // TODO discuss why ok() here ? pretty sure it has no effect)
                writeln!(writer, "{}", new_header).ok();
                writeln!(writer, "{}", seq_str).ok();

            } else {
                //if not found in the map, writing the original
                writeln!(writer, "{}", full_header).ok();
                writeln!(writer, "{}", seq_str).ok();
            }
        }
    })?;

    writer.flush()?; 
    Ok(())
}



// === MERGE ===

// --- MAIN FUNCTION ---

// merge an arbitrary number of indexes  listed in a fof
// each line of the fof is expected to he path to one index dir
// in the end, alsos write a  new metadata CSV file in the output dir
// pub fn merge_multiple_indexes(indexes_fof: &str, output_dir: Option<&str>) -> io::Result<()> {
//     // read the list of index directories.
//     let index_dirs: Vec<String> = {
//         let file = File::open(indexes_fof)?;
//         io::BufReader::new(file)
//             .lines()
//             .filter_map(Result::ok)
//             .map(|line| line.trim().to_string())
//             .filter(|line| !line.is_empty())
//             .collect()
//     };

//     if index_dirs.is_empty() {
//         return Err(io::Error::new(
//             io::ErrorKind::InvalidInput,
//             "No index directories provided in the file-of-indexes",
//         ));
//     }

//     // read metadata from the first index as base parameter
//     let (k, m, bf_size, partition_number, first_color, abundance_number, abundance_max, dense_option) =
//     read_partition_from_csv(&index_dirs[0], "index_info.csv")?;

//     //  vector to store (index_dir, color_count) for every index.
//     let mut indexes_metadata = vec![(index_dirs[0].clone(), first_color)];
//     let mut new_color_number = first_color;

//     // for all other indexes, check that the parameters match and add its color count
//     for index_dir in index_dirs.iter().skip(1) {
//         let (k2, m2, bf_size2, partition_number2, color_count, abundance_number2, abundance_max2, dense_option) =
//         read_partition_from_csv(index_dir, "index_info.csv")?;
//         if k != k2 || m != m2 || bf_size != bf_size2 || partition_number != partition_number2 ||
//            abundance_number != abundance_number2 || abundance_max != abundance_max2 {
//             return Err(io::Error::new(
//                 io::ErrorKind::InvalidInput,
//                 format!("Index {} does not match parameters of the first index", index_dir),
//             ));
//         }
//         new_color_number += color_count;
//         indexes_metadata.push((index_dir.clone(), color_count));
//     }

//     let out_dir = output_dir.unwrap_or("merged_index");
//     fs::create_dir_all(out_dir)?;

//     let partitioned_bf_size = (bf_size as usize) / partition_number;

//     // for each partition, merge the corresponding bfs for every index
//     for partition_idx in 0..partition_number {
//         // for each index, build the file name for the given partition
//         let chunk_files: Vec<String> = indexes_metadata
//             .iter()
//             .map(|(index_dir, _)| {
//                 format!("{}/partition_bloom_filters_p{}.txt", index_dir, partition_idx)
//             })
//             .collect();

//         // collect the color counts from each index
//         let color_counts: Vec<usize> = indexes_metadata.iter().map(|(_, count)| *count).collect();
//         let total_nb_colors = new_color_number;

//         //merge
//         let mut merged_bf = RoaringBitmap::new();
//         merge_partition_bloom_filters(
//             chunk_files,
//             partition_idx,
//             partitioned_bf_size,
//             abundance_number,
//             &color_counts,
//             &mut merged_bf,
//             out_dir,
//             total_nb_colors,
//         )?;
//     }

//     // write the new merged metadata
//     let csv_path = format!("{}/index_info.csv", out_dir);
//     let mut csv_writer = Writer::from_writer(BufWriter::new(File::create(&csv_path)?));
//     csv_writer.write_record(&[
//         "k", "m", "bf_size", "partition_number", "color_number", "abundance_number", "abundance_max"
//     ])?;
//     csv_writer.write_record(&[
//         k.to_string(),
//         m.to_string(),
//         bf_size.to_string(),
//         partition_number.to_string(),
//         new_color_number.to_string(),
//         abundance_number.to_string(),
//         abundance_max.to_string(),
//     ])?;
//     csv_writer.flush()?;

//     println!("Successfully merged {} indexes into directory: {}", indexes_metadata.len(), out_dir);
//     Ok(())
// }

// // --- MERGE FUNCTIONS ---

fn merge_all_partitions(
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
    (0..num_partitions).into_par_iter().try_for_each(|partition_idx| {

        // collect chunk files for the current partition
        let chunk_files_for_partition: Vec<String> = color_counts_per_chunk
            .iter()
            .enumerate()
            .map(|(chunk_idx, _)| {
                format!(
                    "{}/partition_bloom_filters_c{}_p{}.txt",
                    chunk_files_dir, chunk_idx, partition_idx
                )
            })
            .collect();

        // bf for this partition
        let mut partition_bf =RoaringBitmap::new();

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
    println!(
        "All partitions merged and written to disk in {:.2?}",
        elapsed_time
    );

    Ok(())
}

fn merge_partition_bloom_filters(
    chunk_files: Vec<String>,
    partition_idx: usize,
    partitioned_bf_size: usize,
    abundance_number: usize,
    color_counts: &Vec<usize>, // number of colors for each chunk
    bloom_filter:  &mut RoaringBitmap,
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
    let output_file_path = format!("{}/partition_bloom_filters_p{}.txt", output_dir, partition_idx);
    let  output_file = File::create(&output_file_path)?;
    let mut writer = BufWriter::new(output_file);
    writer.write_all(&total_nb_colors.to_le_bytes())?;
    final_bf.serialize_into(&mut writer)?;
    // update the bloom filter for the partition
    *bloom_filter = final_bf;
    // Write the partition's Bloom filter to disk
    

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
                    eprintln!("Failed to load Bloom filter {}: {}", chunk_file, e);
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

                // collect positions in the slice and adjust by offset
                // final_positions.extend(
                //     chunk_bf.clone()
                //         .into_range(slice_start as u32..slice_end as u32)
                //         .map(|pos| pos - slice_start as u32 + current_offset as u32),
                // );
                let slice_start_u32 = slice_start as u32;
                let current_offset_u32 = current_offset as u32;
                final_positions.extend(
                    chunk_bf
                        .range(slice_start_u32..slice_end as u32)
                        .map(|pos| pos - slice_start_u32 + current_offset_u32),
                );

                // update the offset for the next chunk in this slice
                current_offset += abundance_number * chunk_color_number;
            } else {
                eprintln!("Skipping chunk {} due to previous load error", chunk_idx);
            }
        }
    }

    RoaringBitmap::from_sorted_iter(final_positions).expect("Attempt to merge with unsorted positions")
}


// === GENERAL === 

// --- BF MANAGEMENT ---

fn compute_base_position(
    kmer_hash: u64,
    partitioned_bf_size: usize,
    color_number: usize,
    abundance_number: usize,
) -> u64 {
    // return the first position of the concerned column in the partition
    let position = kmer_hash % (partitioned_bf_size as u64);
    position * (color_number as u64) * (abundance_number as u64)
}

// TOUN
/// update color abundances for a specific base position in the Bloom filter
fn update_color_abundances(
    bitmap: &RoaringBitmap,
    base_position: u64,
    color_number: usize,
    abundance_number: usize,
    color_abundances: &mut Vec<Vec<usize>>,
) {
    for color in 0..color_number {
        let mut insert=false;
        for abundance in 0..abundance_number {
            let position_to_check = base_position
                + (color as u64) * (abundance_number as u64)
                + (abundance as u64);
            
            if bitmap.contains(position_to_check as u32) {
                color_abundances[color].push(abundance);
                insert=true;
                break; // keep the minimum
            }
        }
        if !insert{
            // TODO weird discuss value
            color_abundances[color].push(666); // important to record absent k-mers, to compute the median value, also, todo test 
        }
    }
}

fn compute_location_filter(
    hash_kmer: u64, 
    partitioned_bf_size: usize, 
    color_number: usize, // the total nb of indexed fastas (in a chunk)
    path_color_number:usize, // the index of the current indexed fasta
    abundance_number: usize, 
    log_abundance: u16,) 
    -> u64 {
    
    // compute the position to write
    (hash_kmer % (partitioned_bf_size as u64)) * (color_number as u64) * (abundance_number as u64) 
        + (path_color_number as u64) * (abundance_number as u64) 
        + (log_abundance as u64)
    /* example
                        c0  c1  c2  c3	      
        color 0   abund 0   0   1   1
                  abund 1   0   0   0
        color 1   abund 0   0   1   0
                  abund 1   0   0   1
        color 2   abund 0   0   0   1
                  abund 1   0   1   0
        with two partitions becomes two files
        f0=010101000000 

        f1=101001100110

        let's say a k-mer x is in color 2 with abund 0, hashed in column 2
                        c0  c1  c2  c3	      
        color 0   abund 0   0   1   1
                  abund 1   0   0   0
        color 1   abund 0   0   1   0
                  abund 1   0   0   1
        color 2   abund 0   0   X   1
                  abund 1   0   1   0
        
        so we must do a modification so f1 becomes
        f1=1=1010X1100110 with X = 1, so index 4 in the vector

        partition number is 1, and the partition size (number of columns) is 2, composed of column 2 and 3. 
        We'll have hash(x)%2=0, which means we're in the first chunk of the vector ((hash_k % partitioned_bf_size)*color_number*abundance_number = 0 here)
        then we want to go up to color 2, so we mus pass through color 0 and 1 that have two abundances each (path_color_number*abundance_number = 2*2 here)
        then go to the right abundance, here 0 that corresponds to log_abundance
        position_to_write = 0 + 2*2 + 0 = index 4 in the vector
    */
}

fn _get_current_chunk_index(i: usize, chunk_sizes: &Vec<usize>,  partition_nb:usize) -> usize {
    let mut cumulative_size = 0;
    
    for (chunk_idx, &_chunk_size) in chunk_sizes.iter().enumerate() {
        cumulative_size += partition_nb;
        if i < cumulative_size {
            return chunk_idx;
        }
    }
    panic!("Index {} out of bounds for chunk sizes {:?}", i, chunk_sizes);
}

fn update_sparse_counts(
    bloom_filters: &Arc<Vec<Mutex<RoaringBitmap>>>,
    atomic_sparse_one_seen: &atomic::AtomicU64,
    atomic_sparse_fp_seen: &atomic::AtomicU64,
    abundance_number: usize,
) -> io::Result<()> {
    let ones_by_partition = Arc::new(Mutex::new(Vec::<(usize,usize)>::new()));
    bloom_filters.par_iter().for_each(|bitmap| {
        let locked_bitmap = bitmap.lock().unwrap();
        let mut color: usize = 0;
        let mut new_color: usize = 0;
        let mut ones: usize = 0;
        let mut fp: usize = 0;
        locked_bitmap.iter().for_each(|value| {
            ones += 1;
            new_color = (value as usize) / abundance_number;
            if new_color != color {
                color = new_color;
            } else {
                fp += 1;
            }
        });
        let mut tmp_vector = ones_by_partition.lock().unwrap();
        tmp_vector.push((ones, fp));
        atomic_sparse_one_seen.fetch_add(ones as u64, atomic::Ordering::Relaxed);
        atomic_sparse_fp_seen.fetch_add(fp as u64, atomic::Ordering::Relaxed);
    });

    let file = File::create("partitions_info.log")?;
    let mut writer = BufWriter::new(file);
    let tmp_vector: std::sync::MutexGuard<'_, Vec<(usize, usize)>> = ones_by_partition.lock().unwrap();
    for (a1, a2) in tmp_vector.iter() {
        writer.write_all(format!("{},{}\n", a1, a2).as_bytes())?;
    }
    writer.flush()?;

    Ok(())
}

fn flush_map_into_bfs (
    partition_kmers: &mut HashMap<usize, Vec<(u64, u16, usize, usize)>>,
    bloom_filters: &Arc<Vec<Mutex<RoaringBitmap>>>,
    partitioned_bf_size: usize,
    color_number: usize,
    abundance_number: usize,
) -> io::Result<()> {
    // this part fills the BFs per paritition
    for (partition_index, kmers) in partition_kmers.drain() {// iterates and empties the hash map when needed
        let mut kmer_hashes_to_update = Vec::new();
        for (kmer_hash, log_abundance, path_num, _chunk_index) in kmers { // select the bit in the BF
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
        bloom_filter.extend(kmer_hashes_to_update);
    };
    Ok(())
}


// --- WRITE INDEX ON DISK ---

fn create_dir_and_files(num_partition: usize, output_dir: &str) -> io::Result<(Vec<String>, String)> {
    println!("Writing partitioned files in directory: {}", output_dir); 
    let output_path = Path::new(output_dir);
    let partition_dir = match output_path.is_relative() {
        true => std::env::current_dir()?.join(output_path),
        false => output_path.to_path_buf(),
    };
    if !(partition_dir.try_exists()?) {
        fs::create_dir_all(&partition_dir)?;
    }
    let file_paths = Vec::with_capacity(num_partition);
    let partition_dir_string = partition_dir.to_string_lossy().into_owned();
    Ok((file_paths, partition_dir_string))
    /* For instance, if this is the complete structure, 
                            c0  c1  c2  c3	      
            color 0   abund 0   0   1   1
                      abund 1   0   0   0
            color 1   abund 0   0   1   0
                      abund 1   0   0   1
            color 2   abund 0   0   0   1
                      abund 1   0   1   0
    and there are 2 partitions, then the first partition contains
    010101000000
    and the second
    101001100110
    the size is nb of colors (3) * nb of abundances (2) * partition size (in nb of columns, here 2) <- vector_size = 12
    so here the fn will write
    000000000000  as a partitioned_bloom_filters
    */
}

fn write_bloom_filters_to_disk(
    dir_path: &str,
    bloom_filters: &Arc<Vec<Mutex<RoaringBitmap>>>,
    chunks: &[usize], //  number of colors for each chunk
    partition_nb: usize,
    chunk_id:usize
) -> io::Result<()> {
        for (i, bitmap) in bloom_filters.iter().enumerate() {
            // let c = get_current_chunk_index(i, chunks, partition_nb); // chunk idx
            let p = i % partition_nb; // TOUN
            let file_path = Path::new(dir_path).join(format!(
                "partition_bloom_filters_c{}_p{}.txt",
                chunk_id, p
            ));
            let file = File::create(&file_path)?;
            let mut writer = BufWriter::new(file);
            let chunk_colors = chunks[chunk_id];
            // write the number of colors as a u64 to the file first
            writer.write_all(&chunk_colors.to_le_bytes())?;
            // serialize the bitmap into the file
            let mut locked_bitmap = bitmap.lock().unwrap();
            locked_bitmap.serialize_into(&mut writer)?;
            locked_bitmap.clear();
            // bitmap.serialize_into(&mut writer)?;
            // bitmap.clear();
        }
    
    Ok(())
}

fn _write_bloom_filters_to_disk_nochunk(
    dir_path: &str,
    bloom_filters: &[RoaringBitmap],
    nb_colors: usize,
) -> io::Result<()> {
    for (i, bitmap) in bloom_filters.iter().enumerate() {
        let file_path = Path::new(dir_path).join(format!("partition_bloom_filters_p{}.txt", i));
        let file = File::create(&file_path)?;
        let mut writer = BufWriter::new(file);
        writer.write_all(&nb_colors.to_le_bytes())?;
        bitmap.serialize_into(&mut writer)?;
    }
    Ok(())
}

fn write_dense_indexes_to_disk(
    dir_path: &str,
    dense_indexes: &Arc<Vec<Mutex<HashMap<u64, Vec<u8>>>>>
 ) -> io::Result<()> {
    for (partition, hashmap) in dense_indexes.iter().enumerate() {
        let file_path = Path::new(dir_path).join(format!("partition_dense_index_p{}.txt", partition));
        let file = File::create(&file_path)?;
        let mut writer = BufWriter::new(file);
        let mut locked_hashmap = hashmap.lock().unwrap();
        let binary_encoded = bincode::serialize(&locked_hashmap.clone()).unwrap();
        writer.write_all(&binary_encoded)?;
        locked_hashmap.clear();
    }
    Ok(())
}

fn write_kmer_counts_to_disk(
    dir_path: &str,
    kmer_counts_vector: &Arc<Mutex<Vec<usize>>>,
 ) -> io::Result<()> {
    let file_path = Path::new(dir_path).join("kmer_counts_per_color.txt");
    let file = File::create(&file_path)?;
    let mut writer = BufWriter::new(file);
    let mut locked_vector = kmer_counts_vector.lock().unwrap();
    let binary_encoded = bincode::serialize(&locked_vector.clone()).unwrap();
    writer.write_all(&binary_encoded)?;
    locked_vector.clear();
    Ok(())
}

fn write_partition_to_csv(
    bf_dir: &str,
    k: usize,
    m: usize,
    bf_size: u64,
    partition_number: usize,
    color_number: usize,
    abundance_number: usize,
    abundance_max: u16,
    dense_option: bool,
) -> io::Result<()> {
    let output_path = format!("{}/index_info.csv", bf_dir);

    let mut csv_writer = Writer::from_writer(BufWriter::new(File::create(&output_path)?));

    csv_writer.write_record(["k", "m", "bf_size", "partition_number", "color_number", "abundance_number", "abundance_max", "dense_option"])?;
    csv_writer.write_record(&[
            k.to_string(),
            m.to_string(),
            bf_size.to_string(),
            partition_number.to_string(),
            color_number.to_string(),
            abundance_number.to_string(),
            abundance_max.to_string(),
            dense_option.to_string(),
        ])?;

    println!("Index information written to {}", output_path);
    Ok(())
}

// --- LOAD FROM DISK ---

/// load a Bloom filter from disk
fn load_bloom_filter(file_path: &str) -> io::Result<(RoaringBitmap, usize)> {
    let mut file = File::open(file_path)?;

    // read the first 8 bytes as a u64 to get the number of colors
    let mut color_buffer = [0u8; 8];
    file.read_exact(&mut color_buffer)?;
    let local_color_nb = u64::from_le_bytes(color_buffer) as usize;

    // Rread the rest of the file to deserialize the Bloom filter
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;
    let bitmap = RoaringBitmap::deserialize_from(&buffer[..])
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to deserialize bitmap"))?;
    Ok((bitmap, local_color_nb))
}

// TODO Box<Vec<u8>> est surprenant
fn load_dense_index(file_path: &str) -> io::Result<HashMap<u64, Box<Vec<u8>>>> {
    let mut file = File::open(file_path)?;

    // Read the rest of the file to deserialize the hashmap
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;
    let hashmap = bincode::deserialize_from(&buffer[..])
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to deserialize hashmap"))?;
    Ok(hashmap)
}
   
/// Rload index metadata from the CSV.
fn read_partition_from_csv(bf_dir: &str, output_csv: &str) -> io::Result<(usize, usize, u64, usize, usize, usize, u16, bool)> {
    let csv_path = format!("{}/{}", bf_dir, output_csv);
    let mut reader = csv::Reader::from_reader(BufReader::new(File::open(csv_path)?));
    let values = reader.records().next().expect("Index CSV is empty")?;
    let k = values[0].parse().unwrap_or(0);
    let m = values[1].parse().unwrap_or(0);
    let bf_size = values[2].parse().unwrap_or(0);
    let partition_number = values[3].parse().unwrap_or(0);
    let color_number = values[4].parse().unwrap_or(0);
    let abundance_number = values[5].parse().unwrap_or(0);
    let abundance_max = values[6].parse().unwrap_or(0);
    let dense_option = values[7].parse().unwrap_or(false);

    Ok((k, m, bf_size, partition_number, color_number, abundance_number, abundance_max, dense_option))
}

fn load_kmer_counts_vector(dir_path: &str) -> io::Result<Vec<usize>> {
    let mut file = File::open(Path::new(dir_path).join("kmer_counts_per_color.txt"))?;

    // Read the rest of the file to deserialize the hashmap
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;
    let counts_vector = bincode::deserialize_from(&buffer[..])
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to deserialize the counts vector"))?;
    Ok(counts_vector)
}


// --- FOF MANAGEMENT ---

pub fn read_fof_file(file_path: &str) -> io::Result<(Vec<String>, usize)> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut file_paths = Vec::new();
    let mut color_number = 0;
    for line in reader.lines() {
        let line = line?;
        file_paths.push(line.trim().to_string());
        color_number += 1;
    }
    Ok((file_paths, color_number))
}

fn split_fof(lines: &[String]) -> io::Result<(Vec<Vec<String>>, Vec<usize>)> {
    let total_colors = lines.len();

    // chunk max size
    let magic_nb_split = 128;
    // Determine split_factor 
    let split_factor = if total_colors < magic_nb_split {
        1
    } else {
        total_colors.div_ceil(magic_nb_split)
    };

    let mut fof_chunks = vec![vec![]; split_factor];
    let mut chunk_sizes = vec![0; split_factor]; // store the number of colors in each chunk

    for (i, line) in lines.iter().enumerate() {
        let chunk_index = i / magic_nb_split; 
        fof_chunks[chunk_index].push(line.clone());
        chunk_sizes[chunk_index] += 1; // increment of file paths
    }

    Ok((fof_chunks, chunk_sizes))
}

// TODO discuss: used to return a Result, but only the OK variant was returnd
pub fn explore_muset_dir(dir_str: &str) -> (String, String, usize) {
    let dir_path = Path::new(dir_str);
    if !dir_path.is_dir() {
        panic!("{} is not a directory", dir_str);
    }
    let unitigs_path = dir_path.join("unitigs.fa");
    if !unitigs_path.exists() {
        panic!("File not found : {:#?}", unitigs_path);
    }
    let abundance_path = dir_path.join("unitigs.abundance.mat");
    if !abundance_path.exists() {
        panic!("File not found : {:#?}", abundance_path);
    }

    let reader = BufReader::new(File::open(&abundance_path).expect(&format!("Failed to open {:#?}", abundance_path)));
    let err = &format!("Failed to read {:#?}", abundance_path);
    let firstline = reader.lines().next().expect(err).expect(err);
    
    // number of color == number of \t in the first line ?
    // TODO check it's the same
    let color_nb = firstline.chars().filter(|c| *c == '\t').count();
    (unitigs_path.to_string_lossy().to_string(), abundance_path.to_string_lossy().to_string(), color_nb)
}


// --- FASTA FILES PARSING ---

// TOUN
// TODO verif lecture fichiers compresses
fn is_gz_file(file_path: &str) -> io::Result<bool> {
    if file_path.ends_with(".gz") {
        return Ok(true);
    }
    // If not by extension, check first two bytes for GZ magic number
    let mut file = File::open(file_path)?;
    let mut magic = [0u8; 2];
    if file.read_exact(&mut magic).is_ok() {
        // TODO magic number
        return Ok(magic == [0x1F, 0x8B]);
    }
    Ok(false)
}

// TOUN
/// detects input zst format by magic bytes
fn is_zst_file(file_path: &str) -> io::Result<bool> {
    let mut file = File::open(file_path)?;
    let mut magic = [0u8; 4];
    file.read_exact(&mut magic)?;
    Ok(magic == [0x28, 0xB5, 0x2F, 0xFD]) 
}

/// reads input decompressing if in zst or gz format
fn read_file(file_path: &str) -> io::Result<Box<dyn BufRead>> {
    /* Marchet C. */
    if is_zst_file(file_path)? {
        let file = File::open(file_path)?;
        let decompressed = decode_all(BufReader::new(file))?;
        Ok(Box::new(BufReader::new(io::Cursor::new(decompressed))))
    } else if is_gz_file(file_path)? {
        let file = File::open(file_path)?;
        let gz = GzDecoder::new(file);
        Ok(Box::new(BufReader::new(gz)))
    } else {
        let file = File::open(file_path)?;
        Ok(Box::new(BufReader::new(file)))
    }
}


// --- ABUNDANCE PARSING ---

enum HeaderType {
    BCalm,
    Logan,
}

fn determine_header_type(header: &str) -> Result<HeaderType, io::Error> {
    if header.contains("km:f:") {
        Ok(HeaderType::BCalm)
    } else if header.contains("ka:f:") {
        Ok(HeaderType::Logan)
    } else {
        Err(io::Error::other(
            "Header does not contain a recognized count field (km:f: or ka:f:)",
        ))
    }
}

fn extract_count(header: &str, header_type: &HeaderType) -> Result<u16, io::Error> {
    match header_type {
        HeaderType::BCalm => extract_count_from_bcalm_header(header),
        HeaderType::Logan => extract_count_from_logan_header(header),
    }
}

fn extract_count_from_bcalm_header(header: &str) -> Result<u16, io::Error> {
    header.split_whitespace()
        .find(|&part| part.starts_with("km:f:"))
        .and_then(|km_part| {
            km_part.trim_start_matches("km:f:").parse::<f32>().ok()
        })
        .map(|float_val| float_val as u16) //any value over max u16 will be clamped to the maximum possible value
        .ok_or_else(|| io::Error::other("km:f: value not found in header or could not be parsed"))
}

fn extract_count_from_logan_header(header: &str) -> Result<u16, io::Error> {
    header.split_whitespace()
        .find(|&part| part.starts_with("ka:f:")) 
        .and_then(|ka_part| {
            ka_part.trim_start_matches("ka:f:").parse::<f32>().ok()
        })
        .map(|float_val| float_val as u16) 
        .ok_or_else(|| io::Error::other("ka:f: value not found in header or could not be parsed"))
}



// --- ABUNDANCE ENCODING ---

fn compute_log_abundance(count_value: u16, base: f64, max: u16) -> u16 {
    let mut count_valuef = count_value as f64;
    let threshold = 1.0 / (base-1.0);
    if count_valuef <= 0.0 || base <= 0.0 {
        panic!("value and base must be greater than 0");
        // count_valuef=1.0;
    }
    if count_value > max {
        count_valuef = max as f64;
    }
    if count_valuef < threshold {
        (count_valuef - 1.0) as u16
    } else {
        (count_valuef.ln() / base.ln() + (base-1.0).ln() / base.ln() + threshold - 1.0) as u16
    }
}

fn approximate_value(log_value: usize, base: f64) -> u16 {
    if base <= 0.0 {
        panic!("base must be greater than 0");
    }
    let threshold = 1.0 / (base-1.0);
    let logf = log_value as f64;
    if logf < threshold {
        (log_value + 1) as u16
    } else {
        (base.powf((logf+1.0)-threshold)*threshold) as u16
    }
}

// TOUN
fn compute_base(abundance_number: usize, abundance_max: u16) -> f64 {
    let abundance_numberf = abundance_number as f64;
    const TOL: f64 = 1e-9;
    let abundance_maxf = abundance_max as f64;

    // TODO pourquoi ce test ? C'est pas possible non ?
    // TODO voir si utiliser assert est possible
    if abundance_numberf <= 0.0 {
        panic!("abundance number must be greater than 0");
    }

    // TODO pourquoi <= et pas == ?
    if abundance_max <= 0 {
        panic!("Maximal abundance must be greater than 0");
    }


    let equation = |b: f64| -> f64{
        if b <= 1.0 + 1.0 / abundance_maxf {
            return f64::INFINITY; // Avoid invalid logarithms
        }
        (abundance_maxf * (b - 1.0)).ln() / b.ln() + 1.0 / (b - 1.0) - abundance_numberf
    };

    // Set search interval for b
    let mut lower_bound = 1.0 + 1.0 / abundance_maxf + TOL; // Slightly above 1 to avoid division by zero
    let mut upper_bound = 100.0; // Arbitrary large number

    while equation(upper_bound) > 0.0 {
        upper_bound *= 2.0; // Expand search space if necessary
    }

    let mut fa = equation(lower_bound);
    while (upper_bound - lower_bound) > 2.0 * TOL {
        let m = (upper_bound + lower_bound) / 2.0;
        let fm = equation(m);
        if fa * fm <= 0.0 {
            upper_bound = m; // On poursuit avec la moitié gauche
        } else {
            lower_bound = m; // On poursuit avec la moitié droite
            fa = fm;
        }
    }
    (lower_bound + upper_bound) / 2.0
}


// --- MINIMIZER ---


////  iterator over (k-mer, minimizer) pairs for a given sequence
pub struct KmerMinimizerIterator<'a> {
    seq: &'a [u8],
    minima: Vec<u64>,  // minimizers per starting position of kmers
    current: usize,   // current k-mer start position
    k: usize,
}

impl<'a> Iterator for KmerMinimizerIterator<'a> {
    type Item = (&'a [u8], u64); // (kmer, iterator)

    fn next(&mut self) -> Option<Self::Item> {
        if self.current <= self.seq.len().saturating_sub(self.k) {  // below or equal to the last valid starting index 
            let kmer = &self.seq[self.current..self.current + self.k];
            let minimizer = self.minima[self.current];
            self.current += 1;
            Some((kmer, minimizer))
        } else {
            None
        }
    }
}


// minimisers things
// returns an iterator over (k-mer, minimizer) pairs from sequence input
fn kmer_minimizers_seq_level<'a>(
    seq: &'a [u8],
    k: usize,
    m: usize,
) -> impl Iterator<Item = (u64, u64)> + 'a {
    assert!(
        seq.len() >= k,
        "Sequence must be at least k bases long (got {} vs k = {})",
        seq.len(),
        k
    );
    assert!(k >= m, "k must be greater than or equal to m");

    //  collects the hash of every m-mer in the sequence w/ rolling hash
    let m_hashes: Vec<u64> = NtHashIterator::new(seq, m)
        .unwrap_or_else(|err| {
            panic!("Error creating NtHashIterator for m-mers: {:?}", err);
        })
        .collect();

    // compute sliding window to find minimums over m-mer hashes
    // for each k-mer starting at position i, the m-mers inside it are:
    //   m_hashes[i .. i + (k - m + 1)]
    let window_size = k - m + 1;
    // number of k-mers is seq.len() - k + 1; we select 1 minimizer per kmer, so there are m_hashes.len() - window_size + 1 minimizers
    let minimizers = m_hashes.len().saturating_sub(window_size) + 1;
    let mut minima = Vec::with_capacity(minimizers);
    let mut deque: VecDeque<usize> = VecDeque::new();

    // process the first window: indices from 0 .. window_size
    for i in 0..window_size { // positions in the window
        while let Some(&back) = deque.back() { // if there's a value at the back of the queue -> back is the last position P recorded in the queue
            if m_hashes[back] > m_hashes[i] { // if the minimizer list at position P is greater than what it is at the current position => remove from the queue
                deque.pop_back(); // remove large values from the back
            } else {
                break;
            }
        }
        deque.push_back(i); // current position i becomes a candidate for being a minimizer
    }
    // ===> the minimum for the first window has its position recorded at the front of the deque


    /* example
    window_size = 5
    m_hashes = [4, 2, 5, 1, 3]       hash values for positions 0 through 4.
    for i = 0
        deque init : [] -> deque.push_back(0) // [0]
    for i = 1
        m_hashes[0] = 4, with m_hashes[1] = 2
        4 > 2, pop index 0 from the deque
        deque.push_back(1) // [1]
    i = 2
        m_hashes[1] = 2 h m_hashes[2] = 5
        2 < 5 => break
        deque.push_back(2) // [1,2]
    i=3
        m_hashes[2] = 5 m_hashes[3] = 1
        5 > 3, pop index 2 from deque// [1]
        m_hashes[1] = 2, m_hashes[3] = 1
        2 > 1, pop index 1 // []
        deque.push_back(3) //[3]
    i=4
        m_hashes[3] = 1 m_hashes[4] = 3
        1 < 3 , break
        deque.push_back(4) //[3,4]
    => minimizer at pos 0 of deque


    */


    if let Some(&front) = deque.front() {
        minima.push(m_hashes[front]);
    }
    // then the idea carries on, we just have to remove the leftmost index that is no longer in the window each time

    //process the rest of the windows
    for i in window_size..m_hashes.len() {
        // remove indices that are now outside the current window
        while let Some(&front) = deque.front() {
            if front <= i - window_size {
                deque.pop_front();
            } else {
                break;
            }
        }
        // remove elements that are larger than the current element
        while let Some(&back) = deque.back() {
            if m_hashes[back] > m_hashes[i] {
                deque.pop_back();
            } else {
                break;
            }
        }
        deque.push_back(i);
        if let Some(&front) = deque.front() {
            minima.push(m_hashes[front]);
        }
    }

    // at this point, the number of minima should equal the number of k-mers !!
    //assert_eq!(minima.len(), seq.len() - k + 1); //todo remove

    let kmer_hash_iter = NtHashIterator::new(seq, k) // hash kmers
    .unwrap_or_else(|err| {
        panic!("Error creating NtHashIterator for k-mers: {:?}", err)
    });
    // return an iterator over (hashed k-mers,corresponding minimizers)
    kmer_hash_iter.zip(minima)
}


// --- MISC ---

fn _display_progress(total_kmers: u64, start_time: Instant) {
    let elapsed_s = start_time.elapsed().as_secs_f64();
    let elapsed_ms = start_time.elapsed().as_millis() as f64;
    let kmers_per_ms = total_kmers as f64 / elapsed_ms;
    println!(
        "Processed: {} k-mers | Time elapsed: {:.2} s | Rate: {:.2} k-mers/ms",
        total_kmers.to_formatted_string(&Locale::en),
        elapsed_s,
        kmers_per_ms
    );
}





























































































































/* TESTS */

#[allow(unused_imports)]
#[cfg(test)]
mod tests {
    
    use super::*;
    use bio::io::fasta;
    use std::io::Cursor;

    #[test]
    fn test_read_fof_file() {    
        let test_file_path = "test_fof.txt";
    
        // Create the test file
        {
            let mut test_file = File::create(test_file_path).expect("Failed to create test FOF file");
            writeln!(test_file, "/home/test/path/file1.fasta").expect("Failed to write to test FOF file");
            writeln!(test_file, "/home/test/path/file2.fasta").expect("Failed to write to test FOF file");
            writeln!(test_file, "file.fasta").expect("Failed to write to test FOF file");
        }
    
        // Call the function being tested
        let (file_paths, _col_nb) = read_fof_file(test_file_path).unwrap();
    
        // Define the expected result
        let expected = vec![
            "/home/test/path/file1.fasta".to_string(),
            "/home/test/path/file2.fasta".to_string(),
            "file.fasta".to_string(),
        ];
    
        // Assertions
        assert_eq!(file_paths, expected);
    
        // Cleanup
        if Path::new(test_file_path).exists() {
            std::fs::remove_file(test_file_path).expect("Failed to remove test FOF file");
        }
    }

    /*
    #[test]
    fn test_is_valid_kmer() {
        assert_eq!(is_valid_kmer(b"ACCCGNTT"), false);
    }
    */

    #[test]
    fn test_extract_count_from_logan_header2() {
        let header = ">0 ka:f:X  L:+:5:+ L:+:5392806:+  L:-:1:-"; 
        let result = extract_count_from_logan_header(header);
        assert!(result.is_err());
    }
    #[test]
    fn test_extract_count_from_logan_header3() {
        let header = ">0 ka:f:12.4   L:+:5:+ L:+:5392806:+  L:-:1:-";
        let expected = 12;
        let result = extract_count_from_logan_header(header).unwrap();
        assert_eq!(result,expected);
    }
    #[test]
    fn test_extract_count_from_logan_header4() {
        let header = ">0 ka:f:65537   L:+:5:+ L:+:5392806:+  L:-:1:-";
        let expected = 65535;
        let result = extract_count_from_logan_header(header).unwrap();
        assert_eq!(result,expected);
    }

    #[test]
    fn test_compute_log_abundance_gives_zero() {
        let result = compute_log_abundance(1, 2.0, 65535);
        assert_eq!(result, 0, "expected log for 1 to be 0");

        let result2 = compute_log_abundance(1, 1.5, 65535);
        assert_eq!(result2, 0, "expected log for 1 to be 0");
    }

    #[test]
    fn test_compute_log_abundance_is_linear() {
        let b = 1.1;
        for n in 1..9 {
            let result = compute_log_abundance(n, b, 65535);
            assert_eq!(result, n-1, "expected abundance function to be linear for the first values");
        }
        assert_eq!(compute_log_abundance(10, b, 65535), compute_log_abundance(11, b, 65535), "expected abundance function not to be linear at the threshold");
    }

    #[test]
    fn test_compute_log_abundance_positive_values() {
        let result = compute_log_abundance(8, 2.0, 65535);
        assert!(((result as f64) - 3.0).abs() < 1e-6, "expected log base 2 of 8 to be 3");
    }
   
    #[test]
    fn test_compute_log_abundance_with_large_values() {
        let result = compute_log_abundance(65535, 2.0, 65535);
        let expected_log = (65535 as f64).log2().floor() as u16; 
        assert_eq!(result, expected_log,"expected log base 2 of 65535 to be {}, but got {}", expected_log, result);
    }

    #[test]
    fn test_compute_log_abundance_above_max() {
        let result = compute_log_abundance(65535, 2.0, 256);
        let expected_log = (256 as f64).log2().floor() as u16;
        assert_eq!(result, expected_log,"expected abundance to be scaled at log base 2 of 256 should have been {}, but got {}", expected_log, result);
    }

    #[test]
    fn test_determine_header_type_and_extract_count() {
        // Test determining header type
        let bcalm_header = ">0 km:f:12.4 L:+:5:+";
        let logan_header = ">0 ka:f:7.8 L:+:5392806:+";
        let unsupported_header = ">0 other:f:10.0 L:+:1:-";
    
        let bcalm_type = determine_header_type(bcalm_header).unwrap();
        assert!(matches!(bcalm_type, HeaderType::BCalm));
    
        let logan_type = determine_header_type(logan_header).unwrap();
        assert!(matches!(logan_type, HeaderType::Logan));
    
        let unsupported_type = determine_header_type(unsupported_header);
        assert!(unsupported_type.is_err());
    
        // Test extracting counts
        let count_from_bcalm = extract_count(bcalm_header, &HeaderType::BCalm).unwrap();
        assert_eq!(count_from_bcalm, 12);
    
        let count_from_logan = extract_count(logan_header, &HeaderType::Logan).unwrap();
        assert_eq!(count_from_logan, 7);
    }

    #[test]
    #[should_panic(expected = "value and base must be greater than 0")]
    fn test_compute_log_abundance_zero_count_value() {
        compute_log_abundance(0, 2.0, 65535);
    }

    #[test]
    #[should_panic(expected = "value and base must be greater than 0")]
    fn test_compute_log_abundance_non_positive_base() {
        compute_log_abundance(8, 0.0, 65535);
    }
    #[test]
    fn test_approximate_value_with_positive_values() {
        let result = approximate_value(3, 2.0);
        assert_eq!(result, 8, "expected 2^3 to be 8");
    }

    #[test]
    #[should_panic(expected = "base must be greater than 0")]
    fn test_approximate_value_with_zero_base() {
        approximate_value(3, 0.0);
    }

    #[test]
    fn test_approximate_value_with_fractional_base() {
        let base = 2.0f64.sqrt();
        // with sqrt(2), when the approximation increase by 2, the abundance should double
        // since log_{sqrt(2)} (x) = 2 log_2(x)
        let result1 = approximate_value(6, base);
        let result2 = approximate_value(8, base);
        let result3 = approximate_value(10, base);
        assert_eq!(result2/2, result1, "check approximation consistency with base sqrt(2)");
        assert_eq!(result3/2, result2, "check approximation consistency with base sqrt(2)");
    }

    #[test]
    fn test_compute_base_with_positive_values() {
        let result = compute_base(16, 1024);
        assert!((result - 1.5635206).abs() < 1e-6, "expected base for 16 with max 1024 to be ~1.56");
    }

    #[test]
    fn test_compute_base_with_large_abundance_number() {
        let result = compute_base(32, 1024);
        assert!((result - 1.218096).abs() < 1e-6, "expected base for 32 with max 1024 to be approximately ~1.22");
    }

    #[test]
    fn test_compute_base_with_small_abundance_number() {
        let result = compute_base(8, 1024);
        assert!((result - 2.740397).abs() < 1e-6, "expected base for 8 with max 1024 to be ~2.74");
    }

    #[test]
    #[should_panic(expected = "abundance number must be greater than 0")]
    fn test_compute_base_with_zero_abundance_number() {
        compute_base(0, 1024);
    }

    #[test]
    #[should_panic(expected = "Maximal abundance must be greater than 0")]
    fn test_compute_base_with_zero_abundance_max() {
        compute_base(16, 0);
    }


    
    #[test]
    fn test_kmer_hash_minimizers() {
        let seq_str = "ACGTACGTACGTACGT";
        let seq_bytes = seq_str.as_bytes();

        let k = 7;
        let m = 3;

        let pairs: Vec<(u64, u64)> = kmer_minimizers_seq_level(seq_bytes, k, m).collect();

        // there should be seq.len() - k + 1 pairs.
        assert_eq!(pairs.len(), seq_bytes.len() - k + 1);

    }
    
  


    #[test]
    fn test_process_fasta_record() {
        let fasta_input = ">0 ka:f:3.3   L:+:5:+ L:+:5392806:+  L:-:1:-\nAGGAGTAGATACCAGAGATAACGATACAGGTGCGA\n";

        let reader = fasta::Reader::new(Cursor::new(fasta_input));
        let result = reader.records().next().expect("failed to read record");

        let base = 2.0; 
        let processed = process_fasta_record(result, base, 65535, &HeaderType::Logan);

        assert!(processed.is_ok(), "processing failed");
        let (seq, log_abundance, _) = processed.unwrap();
        let expected_seq = b"AGGAGTAGATACCAGAGATAACGATACAGGTGCGA".to_vec();
        let expected_log_abundance = 1; // log2(3) rounds to 1

        assert_eq!(seq, expected_seq, "sequence mismatch");
        assert_eq!(log_abundance, expected_log_abundance, "log abundance mismatch");
    }
   
    
/*
    #[test]
    fn test_update_bloom_filter_memory() {
        // Wrap RoaringBitmap in a Mutex
        let bloom_filter = std::sync::Mutex::new(RoaringBitmap::new());
    
        let kmer_hashes = vec![42, 84, 126]; 
        let partitioned_bf_size = 16;
        let color_number = 3; 
        let path_num = 1;
        let abundance_number = 2; 
        let log_abundance = 1; 
    
        let expected_positions: Vec<u32> = kmer_hashes
            .iter()
            .map(|&hash| {
                compute_location_filter(
                    hash,
                    partitioned_bf_size,
                    color_number,
                    path_num,
                    abundance_number,
                    log_abundance,
                ) as u32
            })
            .collect();
    
        // Pass a reference to the Mutex<RoaringBitmap> instead of &mut RoaringBitmap
        update_bloom_filter_memory(
            &bloom_filter,
            &kmer_hashes,
            partitioned_bf_size,
            color_number,
            path_num,
            abundance_number,
            log_abundance,
        );
    
        // Lock the mutex before checking the bitmap
        let bloom_filter_guard = bloom_filter.lock().unwrap();
    
        for &expected_position in &expected_positions {
            assert!(
                bloom_filter_guard.contains(expected_position),
                "Expected position {} not found in bloom filter",
                expected_position
            );
        }
    
        assert_eq!(
            bloom_filter_guard.len(),
            expected_positions.len().try_into().unwrap(),
            "Unexpected number of entries in bloom filter"
        );
    }
*/
    #[test]
    fn test_update_color_abundances() {
        use roaring::RoaringBitmap;

        let mut bitmap = RoaringBitmap::new();
        let base_position = 100;
        let color_number = 3; 
        let abundance_number = 2; 

        bitmap.insert((base_position + 0) as u32); // color 0, abundance 0
        bitmap.insert((base_position + 1) as u32); // color 0, abundance 1
        bitmap.insert((base_position + 2) as u32); // color 1, abundance 0

        let mut color_abundances = vec![vec![]; color_number];

        update_color_abundances(&bitmap, base_position, color_number, abundance_number, &mut color_abundances);

        let expected_color_abundances = vec![
            vec![0], // color 0 has abundance levels 0 and 1 -> will keep the min
            vec![0],    // color 1 has abundance level 0
            vec![666],     // color 2 has no abundance, set to 666 instead (handle later in the pipeline)
        ];

        assert_eq!(
            color_abundances, expected_color_abundances,
            "Color abundances mismatch: expected {:?}, got {:?}",
            expected_color_abundances, color_abundances
        );
    }


    #[test]
    fn test_write_and_read_partition_csv() {
        let bf_dir = "test_bf_dir";
        let output_csv = "index_info.csv";
        let k = 31;
        let m = 15;
        let bf_size = 1024;
        let partition_number = 4;
        let color_number = 3;
        let abundance_number = 2;
        let abundance_max = 512;
        let dense_option = false;

        fs::create_dir_all(bf_dir).expect("Failed to create test directory");

        let write_result = write_partition_to_csv(
            bf_dir,
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            dense_option,
        );
        assert!(write_result.is_ok(), "Failed to write CSV");

        let read_result = read_partition_from_csv(bf_dir, output_csv);
        assert!(read_result.is_ok(), "Failed to read CSV");

        let (read_k, read_m, read_bf_size, read_partition_number, read_color_number, read_abundance_number, read_abundance_max, read_dense_option) =
            read_result.unwrap();

        assert_eq!(read_k, k, "Mismatch in k value");
        assert_eq!(read_m, m, "Mismatch in m value");
        assert_eq!(read_bf_size, bf_size, "Mismatch in bf_size value");
        assert_eq!(read_partition_number, partition_number, "Mismatch in partition_number value");
        assert_eq!(read_color_number, color_number, "Mismatch in color_number value");
        assert_eq!(read_abundance_number, abundance_number, "Mismatch in abundance_number value");
        assert_eq!(read_abundance_max, abundance_max, "Mismatch in abundance_max value");
        assert_eq!(read_dense_option, dense_option, "Mismatch in dense_option value");

        fs::remove_dir_all(bf_dir).expect("Failed to remove test directory");
    }


    #[test]
    fn test_build_and_query_index_single_sparse() {
        let test_dir = "test_files_bq0";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");
    
        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);
    
        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
        }
        
        {
            let mut file1 = File::create(&file1_path).expect("Failed to create file1.fasta");
            writeln!(file1, ">seq1 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACACAGATCA").expect("Failed to write sequence");
            writeln!(file1, ">seq2 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT").expect("Failed to write sequence");
            writeln!(file1, ">seq3 ka:f:2").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA").expect("Failed to write sequence");
        }
    
        let k = 31;
        let m = 15;
        let bf_size = 1024 * 1024;
        let partition_number = 4;
        let color_number = 1;
        let abundance_number = 256; 
        let abundance_max = 65535;
        let dense_option = false;
        let threshold = color_number;
    
        let (_file_paths, index_dir) = build_index(
            vec![file1_path.clone()],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");
    
        let query_results_path = format!("{}/query_results.csv", index_dir);
        
        query_index(&file1_path, &index_dir, &query_results_path, false, false, 0.5)
            .expect("Failed to query sequences");
    
        let mut reader = csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));
        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let header = record[0].to_string(); // Now reading the header instead of the sequence ID
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((header, color, abundance));
        }
        
        assert!(
            !results.is_empty(),
            "Empty results"
        );
    
        let mut expected_results = vec![
            (">seq1 ka:f:30".to_string(), 0, 29),
            (">seq2 ka:f:30".to_string(), 0, 29),
            (">seq3 ka:f:2".to_string(), 0, 2),  // Values with errors due to log conversion
        ];
    
        results.sort();
        expected_results.sort();



        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }
    
        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");
    }

    #[test]
    fn test_build_and_query_index_single_dense() {
        let test_dir = "test_files_bq0d";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");
    
        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);

    
        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
        }
        
        {
            let mut file1 = File::create(&file1_path).expect("Failed to create file1.fasta");
            writeln!(file1, ">seq1 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACACAGATCA").expect("Failed to write sequence");
            writeln!(file1, ">seq2 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT").expect("Failed to write sequence");
            writeln!(file1, ">seq3 ka:f:2").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA").expect("Failed to write sequence");
        }
    
        let k = 31;
        let m = 15;
        let bf_size = 1024 * 1024;
        let partition_number = 4;
        let color_number = 1;
        let abundance_number = 255; 
        let abundance_max = 65535;
        let dense_option = true;
        let threshold = color_number;
    
        let (_file_paths, index_dir) = build_index(
            vec![file1_path.clone()],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");
    
        let query_results_path = format!("{}/query_results.csv", index_dir);
        
        query_index(&file1_path, &index_dir, &query_results_path, false, false, 0.5)
            .expect("Failed to query sequences");
    
        let mut reader = csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));
        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let header = record[0].to_string(); // Now reading the header instead of the sequence ID
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((header, color, abundance));
        }
        
        assert!(
            !results.is_empty(),
            "Empty results"
        );
    
        let mut expected_results = vec![
            (">seq1 ka:f:30".to_string(), 0, 29),
            (">seq2 ka:f:30".to_string(), 0, 29),
            (">seq3 ka:f:2".to_string(), 0, 2),  // Values with errors due to log conversion
        ];
    
        results.sort();
        expected_results.sort();



        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }
    
        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");
    }

    #[test]
    fn test_build_and_query_index_sparse() {
        let test_dir = "test_files_bq1";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");
    
        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);
        let file2_path = format!("{}/file2Q.fa", test_dir);
    
        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file2_path).expect("Failed to write to fof.txt");
        }
        
        {
            let mut file1 = File::create(&file1_path).expect("Failed to create file1.fasta");
            writeln!(file1, ">seq1 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACACAGATCA").expect("Failed to write sequence");
            writeln!(file1, ">seq2 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT").expect("Failed to write sequence");
            writeln!(file1, ">seq3 ka:f:2").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA").expect("Failed to write sequence");
        }
    
        {
            let mut file2 = File::create(&file2_path).expect("Failed to create file2.fasta");
            writeln!(file2, ">seq4 ka:f:1000").expect("Failed to write header");
            writeln!(file2, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG").expect("Failed to write sequence");
            writeln!(file2, ">seq5 ka:f:1000").expect("Failed to write header");
            writeln!(file2, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA").expect("Failed to write sequence");
        }
        
        let k = 31;
        let m = 15;
        let bf_size = 1024 * 1024;
        let partition_number = 4;
        let color_number = 2;
        let abundance_number = 256; 
        let abundance_max = 65535;
        let dense_option = false;
        let threshold = color_number;
    
        let (_file_paths, index_dir) = build_index(
            vec![file1_path.clone(), file2_path.clone()],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");
    
        let query_results_path = format!("{}/query_results.csv", index_dir);
        
        query_index(&file1_path, &index_dir, &query_results_path, false, false, 0.5)
            .expect("Failed to query sequences");
    
        let mut reader = csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));
        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let header = record[0].to_string(); // Now reading the header instead of the sequence ID
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((header, color, abundance));
        }
        
        assert!(
            !results.is_empty(),
            "Empty results"
        );
    
        let mut expected_results = vec![
            (">seq1 ka:f:30".to_string(), 0, 29),
            (">seq2 ka:f:30".to_string(), 0, 29),
            (">seq3 ka:f:2".to_string(), 0, 2), 
            (">seq3 ka:f:2".to_string(), 1, 997),
        ];
        results.sort();
        expected_results.sort();
        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }
    
        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");
    }

    #[test]
    fn test_build_and_query_index_dense() {
        let test_dir = "test_files_bq1d";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");
    
        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);
        let file2_path = format!("{}/file2Q.fa", test_dir);
    
        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file2_path).expect("Failed to write to fof.txt");
        }
        
        {
            let mut file1 = File::create(&file1_path).expect("Failed to create file1.fasta");
            writeln!(file1, ">seq1 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACACAGATCA").expect("Failed to write sequence");
            writeln!(file1, ">seq2 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT").expect("Failed to write sequence");
            writeln!(file1, ">seq3 ka:f:2").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA").expect("Failed to write sequence");
        }
    
        {
            let mut file2 = File::create(&file2_path).expect("Failed to create file2.fasta");
            writeln!(file2, ">seq4 ka:f:1000").expect("Failed to write header");
            writeln!(file2, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG").expect("Failed to write sequence");
            writeln!(file2, ">seq5 ka:f:1000").expect("Failed to write header");
            writeln!(file2, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA").expect("Failed to write sequence");
        }
        
        let k = 31;
        let m = 15;
        let bf_size = 1024 * 1024;
        let partition_number = 4;
        let color_number = 2;
        let abundance_number = 256; 
        let abundance_max = 65535;
        let dense_option = true;
        let threshold = color_number;
    
        let (_file_paths, index_dir) = build_index(
            vec![file1_path.clone(), file2_path.clone()],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");
    
        let query_results_path = format!("{}/query_results.csv", index_dir);
        
        query_index(&file1_path, &index_dir, &query_results_path, false, false, 0.5)
            .expect("Failed to query sequences");
    
        let mut reader = csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));
        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let header = record[0].to_string(); // Now reading the header instead of the sequence ID
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((header, color, abundance));
        }
        
        assert!(
            !results.is_empty(),
            "Empty results"
        );
    
        let mut expected_results = vec![
            (">seq1 ka:f:30".to_string(), 0, 29),
            (">seq2 ka:f:30".to_string(), 0, 29),
            (">seq3 ka:f:2".to_string(), 0, 2), 
            (">seq3 ka:f:2".to_string(), 1, 997),
        ];
        results.sort();
        expected_results.sort();
        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }
    
        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");
    }
    

    #[test]
    fn test_build_and_query_index_longseq_simple_sparse() {
        let test_dir = "test_files_bq_ls";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);
        let file2_path = format!("{}/file2Q.fa", test_dir);

        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file2_path).expect("Failed to write to fof.txt");
        }

        {
            let mut file1 = File::create(&file1_path).expect("Failed to create file1.fasta");
            writeln!(file1, ">seq1 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACACAGATCA").expect("Failed to write sequence");
            writeln!(file1, ">seq2 ka:f:10").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT").expect("Failed to write sequence");
        }

        {
            let mut file2 = File::create(&file2_path).expect("Failed to create file2.fasta");
            writeln!(file2, ">seq3 ka:f:1000").expect("Failed to write header");
            writeln!(
                file2,
                "AAAAATGATAGTAGAAAAAAATTTTAAAAAAACACCCCTGG"
            )
            .expect("Failed to write sequence");
        }

        let k = 31;
        let m = 15;
        let bf_size = 1024 * 1024;
        let partition_number = 4;
        let color_number = 2;
        let abundance_number = 256;
        let abundance_max = 65535;
        let dense_option = false;
        let threshold = color_number;

        let (_file_paths, index_dir) = build_index(
            vec![file1_path.clone(), file2_path.clone()],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");

        let query_results_path = format!("{}/query_results.csv", index_dir);

        query_index(&file2_path, &index_dir, &query_results_path, false, false, 0.5)
            .expect("Failed to query sequences");

        // Validate the results written to the query results CSV file
        let mut reader =
            csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));
        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let header = record[0].to_string(); // Now reading the header instead of the sequence ID
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((header, color, abundance));
        }

        assert!(!results.is_empty(), "Empty results");

        let mut expected_results = vec![(">seq3 ka:f:1000".to_string(), 1, 997)];
        results.sort();
        expected_results.sort();
        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }

        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");
    }

    #[test]
    fn test_build_and_query_index_longseq_simple_dense() {
        let test_dir = "test_files_bq_lsd";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);
        let file2_path = format!("{}/file2Q.fa", test_dir);

        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file2_path).expect("Failed to write to fof.txt");
        }

        {
            let mut file1 = File::create(&file1_path).expect("Failed to create file1.fasta");
            writeln!(file1, ">seq1 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACACAGATCA").expect("Failed to write sequence");
            writeln!(file1, ">seq2 ka:f:10").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT").expect("Failed to write sequence");
        }

        {
            let mut file2 = File::create(&file2_path).expect("Failed to create file2.fasta");
            writeln!(file2, ">seq3 ka:f:1000").expect("Failed to write header");
            writeln!(
                file2,
                "AAAAATGATAGTAGAAAAAAATTTTAAAAAAACACCCCTGG"
            )
            .expect("Failed to write sequence");
        }

        let k = 31;
        let m = 15;
        let bf_size = 1024 * 1024;
        let partition_number = 4;
        let color_number = 2;
        let abundance_number = 256;
        let abundance_max = 65535;
        let dense_option = true;
        let threshold = color_number;

        let (_file_paths, index_dir) = build_index(
            vec![file1_path.clone(), file2_path.clone()],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");

        let query_results_path = format!("{}/query_results.csv", index_dir);

        query_index(&file2_path, &index_dir, &query_results_path, false, false, 0.5)
            .expect("Failed to query sequences");

        // Validate the results written to the query results CSV file
        let mut reader =
            csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));
        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let header = record[0].to_string(); // Now reading the header instead of the sequence ID
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((header, color, abundance));
        }

        assert!(!results.is_empty(), "Empty results");

        let mut expected_results = vec![(">seq3 ka:f:1000".to_string(), 1, 997)];
        results.sort();
        expected_results.sort();
        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }

        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");
    }


    #[test]
    fn test_build_and_query_index2_sparse() {
        use std::io::Write;

        let test_dir = "test_files_bq2";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);
        let file2_path = format!("{}/file2Q.fa", test_dir);

        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file2_path).expect("Failed to write to fof.txt");
        }
    
        {
            let mut file1 = File::create(&file1_path).expect("Failed to create file1.fasta");
            writeln!(file1, ">seq1 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACACAGATCA").expect("Failed to write sequence");
            writeln!(file1, ">seq2 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT").expect("Failed to write sequence");
            writeln!(file1, ">seq3 ka:f:2").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA").expect("Failed to write sequence");
        }

        {
            let mut file2 = File::create(&file2_path).expect("Failed to create file2.fasta");
            writeln!(file2, ">seq4 ka:f:1000").expect("Failed to write header");
            writeln!(file2, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG").expect("Failed to write sequence");
            writeln!(file2, ">seq5 ka:f:1000").expect("Failed to write header");
            writeln!(file2, "AAAAAAAAAAAAAAAAAAAAACACCCCTGGG").expect("Failed to write sequence");
        }
     
        let k = 31;
        let m = 15;
        let bf_size = 1024 * 1024;
        let partition_number = 4;
        let color_number = 2;
        let abundance_number = 255; 
        let abundance_max = 65535;
        let dense_option = false;
        let threshold = 1;

       

        let (_file_paths, index_dir) = build_index(
            vec![file1_path.clone(), file2_path.clone()],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");
     
        /*
        fs::rename("test_files_bq2/partitioned_bloom_filters_c0_p0.txt", "test_files_bq2/partitioned_bloom_filters_p0.txt");
        fs::rename("test_files_bq2/partitioned_bloom_filters_c0_p1.txt", "test_files_bq2/partitioned_bloom_filters_p1.txt");
        fs::rename("test_files_bq2/partitioned_bloom_filters_c0_p2.txt", "test_files_bq2/partitioned_bloom_filters_p2.txt");
        fs::rename("test_files_bq2/partitioned_bloom_filters_c0_p3.txt", "test_files_bq2/partitioned_bloom_filters_p3.txt");
        */
        let query_results_path = format!("{}/query_results.csv", index_dir);
        query_index(&file2_path, test_dir, &query_results_path, false, false, 0.5) // query sequences from file1
            .expect("Failed to query sequences");

        
        let mut reader = csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));

        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let seq_id = record[0].to_string();
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((seq_id, color, abundance));
        }

        assert!(
            !results.is_empty(),
            "Empty results"
        );

        let mut expected_results = vec![
            (">seq4 ka:f:1000".to_string(), 1, 979),
            (">seq5 ka:f:1000".to_string(), 1, 979),
        ];
        results.sort();
        expected_results.sort();
        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }
        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");
    }

    #[test]
    fn test_build_and_query_index2_dense() {
        use std::io::Write;

        let test_dir = "test_files_bq2d";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);
        let file2_path = format!("{}/file2Q.fa", test_dir);

        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file2_path).expect("Failed to write to fof.txt");
        }
    
        {
            let mut file1 = File::create(&file1_path).expect("Failed to create file1.fasta");
            writeln!(file1, ">seq1 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACACAGATCA").expect("Failed to write sequence");
            writeln!(file1, ">seq2 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT").expect("Failed to write sequence");
            writeln!(file1, ">seq3 ka:f:2").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA").expect("Failed to write sequence");
        }

        {
            let mut file2 = File::create(&file2_path).expect("Failed to create file2.fasta");
            writeln!(file2, ">seq4 ka:f:1000").expect("Failed to write header");
            writeln!(file2, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG").expect("Failed to write sequence");
            writeln!(file2, ">seq5 ka:f:1000").expect("Failed to write header");
            writeln!(file2, "AAAAAAAAAAAAAAAAAAAAACACCCCTGGG").expect("Failed to write sequence");
        }
     
        let k = 31;
        let m = 15;
        let bf_size = 1024 * 1024;
        let partition_number = 4;
        let color_number = 2;
        let abundance_number = 255; 
        let abundance_max = 65535;
        let dense_option = true;
        let threshold = 1;

       

        let (_file_paths, index_dir) = build_index(
            vec![file1_path.clone(), file2_path.clone()],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");
     
        /*
        fs::rename("test_files_bq2/partitioned_bloom_filters_c0_p0.txt", "test_files_bq2/partitioned_bloom_filters_p0.txt");
        fs::rename("test_files_bq2/partitioned_bloom_filters_c0_p1.txt", "test_files_bq2/partitioned_bloom_filters_p1.txt");
        fs::rename("test_files_bq2/partitioned_bloom_filters_c0_p2.txt", "test_files_bq2/partitioned_bloom_filters_p2.txt");
        fs::rename("test_files_bq2/partitioned_bloom_filters_c0_p3.txt", "test_files_bq2/partitioned_bloom_filters_p3.txt");
        */
        let query_results_path = format!("{}/query_results.csv", index_dir);
        query_index(&file2_path, test_dir, &query_results_path, false, false, 0.5) // query sequences from file1
            .expect("Failed to query sequences");

        
        let mut reader = csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));

        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let seq_id = record[0].to_string();
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((seq_id, color, abundance));
        }

        assert!(
            !results.is_empty(),
            "Empty results"
        );

        let mut expected_results = vec![
            (">seq4 ka:f:1000".to_string(), 1, 979),
            (">seq5 ka:f:1000".to_string(), 1, 979),
        ];
        results.sort();
        expected_results.sort();
        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }
        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");
    }


    #[test]
    fn test_build_and_query_index_sharedk_sparse() {

        let test_dir = "test_files_bq3";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);
        let file2_path = format!("{}/file2Q.fa", test_dir);
        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file2_path).expect("Failed to write to fof.txt");
        }

        {
            let mut file1 = File::create(&file1_path).expect("Failed to create file1.fasta");
            writeln!(file1, ">seq1 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACACAGATCA").expect("Failed to write sequence");
            writeln!(file1, ">seq2 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT").expect("Failed to write sequence");
            writeln!(file1, ">seq3 ka:f:1500").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA").expect("Failed to write sequence");
        }

        {
            let mut file2 = File::create(&file2_path).expect("Failed to create file2.fasta");
            writeln!(file2, ">seq4 ka:f:1000").expect("Failed to write header");
            writeln!(file2, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG").expect("Failed to write sequence");
            writeln!(file2, ">seq5 ka:f:1000").expect("Failed to write header");
            writeln!(file2, "AAAAAAAAAAAAAAAAAAAAACACCCCTGGG").expect("Failed to write sequence");
            writeln!(file2, ">seq6 ka:f:4").expect("Failed to write header");
            writeln!(file2, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA").expect("Failed to write sequence");
        }
        
        let k = 31;
        let m = 15;
        let bf_size = 1024 * 1024;
        let partition_number = 4;
        let color_number = 2;
        let abundance_number = 255; 
        let abundance_max = 255;
        let dense_option = false;
        let threshold = color_number;


        let (_file_paths, index_dir) = build_index(
            vec![file1_path.clone(), file2_path.clone()],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");
        /*
        fs::rename("test_files_bq3/partitioned_bloom_filters_c0_p0.txt", "test_files_bq3/partitioned_bloom_filters_p0.txt");
        fs::rename("test_files_bq3/partitioned_bloom_filters_c0_p1.txt", "test_files_bq3/partitioned_bloom_filters_p1.txt");
        fs::rename("test_files_bq3/partitioned_bloom_filters_c0_p2.txt", "test_files_bq3/partitioned_bloom_filters_p2.txt");
        fs::rename("test_files_bq3/partitioned_bloom_filters_c0_p3.txt", "test_files_bq3/partitioned_bloom_filters_p3.txt");
        */
        let query_results_path = format!("{}/query_results.csv", index_dir);
        query_index(&file1_path, test_dir, &query_results_path, false, false, 0.5) // query sequences from file1
            .expect("Failed to query sequences");

        
        let mut reader = csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));

        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let seq_id = record[0].to_string();
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((seq_id, color, abundance));
        }

        assert!(
            !results.is_empty(),
            "Empty results"
        );
        
        let mut expected_results = vec![
            (">seq1 ka:f:30".to_string(), 0, 30),
            (">seq2 ka:f:30".to_string(), 0, 30),
            (">seq3 ka:f:1500".to_string(), 0, 255), //because of abundance_max
            (">seq3 ka:f:1500".to_string(), 1, 4),
        ];
        results.sort();
        expected_results.sort();
        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }
        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");
            
    }

    #[test]
    fn test_build_and_query_index_sharedk_dense() {

        let test_dir = "test_files_bq3d";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);
        let file2_path = format!("{}/file2Q.fa", test_dir);
        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file2_path).expect("Failed to write to fof.txt");
        }

        {
            let mut file1 = File::create(&file1_path).expect("Failed to create file1.fasta");
            writeln!(file1, ">seq1 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACACAGATCA").expect("Failed to write sequence");
            writeln!(file1, ">seq2 ka:f:30").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT").expect("Failed to write sequence");
            writeln!(file1, ">seq3 ka:f:1500").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA").expect("Failed to write sequence");
        }

        {
            let mut file2 = File::create(&file2_path).expect("Failed to create file2.fasta");
            writeln!(file2, ">seq4 ka:f:1000").expect("Failed to write header");
            writeln!(file2, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG").expect("Failed to write sequence");
            writeln!(file2, ">seq5 ka:f:1000").expect("Failed to write header");
            writeln!(file2, "AAAAAAAAAAAAAAAAAAAAACACCCCTGGG").expect("Failed to write sequence");
            writeln!(file2, ">seq6 ka:f:4").expect("Failed to write header");
            writeln!(file2, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA").expect("Failed to write sequence");
        }
        
        let k = 31;
        let m = 15;
        let bf_size = 1024 * 1024;
        let partition_number = 4;
        let color_number = 2;
        let abundance_number = 255; 
        let abundance_max = 255;
        let dense_option = true;
        let threshold = color_number;


        let (_file_paths, index_dir) = build_index(
            vec![file1_path.clone(), file2_path.clone()],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");
        /*
        fs::rename("test_files_bq3/partitioned_bloom_filters_c0_p0.txt", "test_files_bq3/partitioned_bloom_filters_p0.txt");
        fs::rename("test_files_bq3/partitioned_bloom_filters_c0_p1.txt", "test_files_bq3/partitioned_bloom_filters_p1.txt");
        fs::rename("test_files_bq3/partitioned_bloom_filters_c0_p2.txt", "test_files_bq3/partitioned_bloom_filters_p2.txt");
        fs::rename("test_files_bq3/partitioned_bloom_filters_c0_p3.txt", "test_files_bq3/partitioned_bloom_filters_p3.txt");
        */
        let query_results_path = format!("{}/query_results.csv", index_dir);
        query_index(&file1_path, test_dir, &query_results_path, false, false, 0.5) // query sequences from file1
            .expect("Failed to query sequences");

        
        let mut reader = csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));

        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let seq_id = record[0].to_string();
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((seq_id, color, abundance));
        }

        assert!(
            !results.is_empty(),
            "Empty results"
        );
        
        let mut expected_results = vec![
            (">seq1 ka:f:30".to_string(), 0, 30),
            (">seq2 ka:f:30".to_string(), 0, 30),
            (">seq3 ka:f:1500".to_string(), 0, 255), //because of abundance_max
            (">seq3 ka:f:1500".to_string(), 1, 4),
        ];
        results.sort();
        expected_results.sort();
        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }
        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");
            
    }

    #[test]
    fn test_build_and_query_index_long_sparse() {

        let test_dir = "test_files_bql";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);
        let file2_path = format!("{}/file2Q.fa", test_dir);
        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file2_path).expect("Failed to write to fof.txt");
        }

        {
            let mut file1 = File::create(&file1_path).expect("Failed to create file1.fasta");
            writeln!(file1, ">seq1 ka:f:30").expect("Failed to write header");
            //writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACACAGATCAGAGGAT").expect("Failed to write sequence");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACACAGATCAGAG").expect("Failed to write sequence");
            writeln!(file1, ">seq2 ka:f:8").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT").expect("Failed to write sequence");
        }

        {
            let mut file2 = File::create(&file2_path).expect("Failed to create file2.fasta");
            writeln!(file2, ">seq4 ka:f:1000").expect("Failed to write header");
            writeln!(file2, "CAAAAAAAAAAAAAAAAAAAAACACCCCTGGAC").expect("Failed to write sequence");
        }
        
        let k = 31;
        let m = 15;
        let bf_size = 1024 * 1024;
        let partition_number = 4;
        let color_number = 2;
        let abundance_number = 255; 
        let abundance_max = 65535;
        let dense_option = false;
        let threshold = color_number;

        let (_file_paths, index_dir) = build_index(
            vec![file1_path.clone(), file2_path.clone()],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");
        /*
        fs::rename("test_files_bq3/partitioned_bloom_filters_c0_p0.txt", "test_files_bq3/partitioned_bloom_filters_p0.txt");
        fs::rename("test_files_bq3/partitioned_bloom_filters_c0_p1.txt", "test_files_bq3/partitioned_bloom_filters_p1.txt");
        fs::rename("test_files_bq3/partitioned_bloom_filters_c0_p2.txt", "test_files_bq3/partitioned_bloom_filters_p2.txt");
        fs::rename("test_files_bq3/partitioned_bloom_filters_c0_p3.txt", "test_files_bq3/partitioned_bloom_filters_p3.txt");
        */
        let query_results_path = format!("{}/query_results.csv", index_dir);
        query_index(&file1_path, test_dir, &query_results_path, false, false, 0.5) // query sequences from file1
            .expect("Failed to query sequences");

       
        let mut reader = csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));

        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let seq_id = record[0].to_string();
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((seq_id, color, abundance));
        }

        assert!(
            !results.is_empty(),
            "Empty results"
        );
        
        let mut expected_results = vec![
            (">seq1 ka:f:30".to_string(), 0, 29),
            (">seq2 ka:f:8".to_string(), 0, 8),
        ];
        results.sort();
        expected_results.sort();
        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }
        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");
            
    }


    #[test]
    fn test_build_and_query_index_long_dense() {

        let test_dir = "test_files_bqld";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);
        let file2_path = format!("{}/file2Q.fa", test_dir);
        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file2_path).expect("Failed to write to fof.txt");
        }

        {
            let mut file1 = File::create(&file1_path).expect("Failed to create file1.fasta");
            writeln!(file1, ">seq1 ka:f:30").expect("Failed to write header");
            //writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACACAGATCAGAGGAT").expect("Failed to write sequence");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACACAGATCAGAG").expect("Failed to write sequence");
            writeln!(file1, ">seq2 ka:f:8").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT").expect("Failed to write sequence");
        }

        {
            let mut file2 = File::create(&file2_path).expect("Failed to create file2.fasta");
            writeln!(file2, ">seq4 ka:f:1000").expect("Failed to write header");
            writeln!(file2, "CAAAAAAAAAAAAAAAAAAAAACACCCCTGGAC").expect("Failed to write sequence");
        }
        
        let k = 31;
        let m = 15;
        let bf_size = 1024 * 1024;
        let partition_number = 4;
        let color_number = 2;
        let abundance_number = 255; 
        let abundance_max = 65535;
        let dense_option = true;
        let threshold = color_number;

        let (_file_paths, index_dir) = build_index(
            vec![file1_path.clone(), file2_path.clone()],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");
        /*
        fs::rename("test_files_bq3/partitioned_bloom_filters_c0_p0.txt", "test_files_bq3/partitioned_bloom_filters_p0.txt");
        fs::rename("test_files_bq3/partitioned_bloom_filters_c0_p1.txt", "test_files_bq3/partitioned_bloom_filters_p1.txt");
        fs::rename("test_files_bq3/partitioned_bloom_filters_c0_p2.txt", "test_files_bq3/partitioned_bloom_filters_p2.txt");
        fs::rename("test_files_bq3/partitioned_bloom_filters_c0_p3.txt", "test_files_bq3/partitioned_bloom_filters_p3.txt");
        */
        let query_results_path = format!("{}/query_results.csv", index_dir);
        query_index(&file1_path, test_dir, &query_results_path, false, false, 0.5) // query sequences from file1
            .expect("Failed to query sequences");

       
        let mut reader = csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));

        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let seq_id = record[0].to_string();
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((seq_id, color, abundance));
        }

        assert!(
            !results.is_empty(),
            "Empty results"
        );
        
        let mut expected_results = vec![
            (">seq1 ka:f:30".to_string(), 0, 29),
            (">seq2 ka:f:8".to_string(), 0, 8),
        ];
        results.sort();
        expected_results.sort();
        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }
        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");

    }

    #[ignore]
    #[test]
    fn test_insert_and_query_kmer_with_verifications() {
        use roaring::RoaringBitmap;
    
        let kmer_hash: u64 = 42; 
        let minimizer: u64 = 7; 
        let partition_number = 4;
        let bf_size = 1024; 
        let partitioned_bf_size = bf_size as usize / partition_number;
        let color_number = 3;
        let path_color_number = 0;
        let abundance_number = 2;
        let log_abundance = 1;
    
        let mut bloom_filters: Vec<RoaringBitmap> = vec![RoaringBitmap::new(); partition_number];
    
        let partition_index_insert = (minimizer % (partition_number as u64)) as usize;
        let position_to_write = compute_location_filter(
            kmer_hash,
            partitioned_bf_size,
            color_number,
            path_color_number,
            abundance_number,
            log_abundance,
        );
        bloom_filters[partition_index_insert].insert(position_to_write as u32);
    
        let nt_hash_iterator = vec![kmer_hash].into_iter();
        let min_iter = vec![(minimizer, 0)].into_iter();
    
        let mut color_abundances = vec![Vec::new(); color_number];
        for (kmer_hash, (minimizer, _position)) in nt_hash_iterator.zip(min_iter) {
            let partition_index_query = (minimizer % (partition_number as u64)) as usize;
    
            assert_eq!(
                partition_index_insert, partition_index_query,
                "Partition index mismatch: insert = {}, query = {}",
                partition_index_insert, partition_index_query
            );
    
            let bitmap = &bloom_filters[partition_index_query];
    
            let base_position = compute_base_position(
                kmer_hash,
                partitioned_bf_size,
                color_number,
                abundance_number,
            );
    
            // I add + log_abund because base_position does not have this info and just finds approximately the value
            assert_eq!(
                position_to_write, base_position+(log_abundance as u64),
                "Position mismatch: insert = {}, query = {}",
                position_to_write, base_position+(log_abundance as u64)
            );
    
            update_color_abundances(
                bitmap,
                base_position,
                color_number,
                abundance_number,
                &mut color_abundances,
            );
        }
    
        let results: Vec<(usize, usize)> = color_abundances
            .into_iter()
            .enumerate()
            .filter_map(|(color, abundances)| {
                let min_abundance = abundances.into_iter().min();
                min_abundance.map(|abundance| (color, abundance))
            })
            .collect();
    
        let expected_results = vec![(0, log_abundance as usize),(1,0),(2,0)]; // expect color 0 with abundance level 1
    
        assert_eq!(
            results, expected_results,
            "Mismatch in query results: expected {:?}, got {:?}",
            expected_results, results
        );
    }



    #[ignore]
    #[test]
    fn test_split_fof_1() -> std::io::Result<()> {
        use std::fs::{self, File};
        use std::io::{BufReader, BufRead, Write};

        let test_dir = "test_files_fofs1";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fof_path = format!("{}/fof_split.txt", test_dir);

        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create test FOF file");
            writeln!(fof_file, "/path/to/file1.fasta").expect("Failed to write to FOF file");
            writeln!(fof_file, "/path/to/file2.fasta").expect("Failed to write to FOF file");
            writeln!(fof_file, "/path/to/file3.fasta").expect("Failed to write to FOF file");
            writeln!(fof_file, "/path/to/file4.fasta").expect("Failed to write to FOF file");
            writeln!(fof_file, "/path/to/file5.fasta").expect("Failed to write to FOF file");
        }

        let file = File::open(&fof_path)?;
        let reader = BufReader::new(file);
        let lines: Vec<String> = reader.lines().collect::<Result<_, _>>()?;

        let result = split_fof(&lines);

        assert!(result.is_ok(), "split_fof returned an error");

        let (fof_chunks, _) = result.unwrap();



        let expected_chunks =
            vec![vec![
                "/path/to/file1.fasta".to_string(),
                "/path/to/file2.fasta".to_string(),]
                ,vec![
                "/path/to/file3.fasta".to_string(),
                "/path/to/file4.fasta".to_string(),]
                ,vec![
                "/path/to/file5.fasta".to_string(),]
            ];

        assert_eq!(
            fof_chunks, expected_chunks,
            "Mismatch in chunk distribution: expected {:?}, got {:?}",
            expected_chunks, fof_chunks
        );

        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");

        Ok(())
        
    }

/*
    #[test]
    fn test_merge_bloom_filters() {
        use roaring::RoaringBitmap;

        // initialize bf1: represents 2 datasets, abundance_number = 3, partitioned_bf_size = 2
        let mut bf1 = RoaringBitmap::new();  //110001,000000
        bf1.insert(0); // 1st dataset, 1st abundance level
        bf1.insert(1); // 1st dataset, 2nd abundance level
        bf1.insert(5); // 2nd dataset, 1st abundance level

        // initialize bf2: represents 1 dataset, abundance_number = 3
        let mut bf2 = RoaringBitmap::new(); // 111,000
        bf2.insert(0); // 1st dataset, 1st abundance level
        bf2.insert(1); // 1st dataset, 2nd abundance level
        bf2.insert(2); // 1st dataset, 3rd abundance level

        let partitioned_bf_size = 2;
        let merged_color_number = 2; // bf1 has 2 colors (datasets)
        let new_color_number = 1;    // bf2 has 1 color (dataset)
        let abundance_number = 3;

        let (merged_bf, color_nb_merge_final) =
            merge_bloom_filters(&bf1, &bf2, partitioned_bf_size, merged_color_number, new_color_number, abundance_number);

        let mut expected_bf = RoaringBitmap::new();
        //  first column: 110001 + 111 -> 110001111
        expected_bf.insert(0);
        expected_bf.insert(1);
        expected_bf.insert(5);
        expected_bf.insert(6);
        expected_bf.insert(7);
        expected_bf.insert(8);
        // second column: 000000 + 000 -> 000000000
        

        let expected_color_nb_merge_final = 3; // 2 (from bf1) + 1 (from bf2)
        assert_eq!(
            color_nb_merge_final, expected_color_nb_merge_final,
            "Final color number does not match the expected result"
        );
        assert_eq!(
            merged_bf, expected_bf,
            "Merged bitmap does not match the expected result"
        );
        
    }
*/


    //#[ignore]
    #[test]
    fn test_build_and_query_index_with_chunks_and_merge() {
        let test_dir = "test_files_bq_merge";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);
        let file2_path = format!("{}/file2Q.fa", test_dir);
        let file3_path = format!("{}/file3Q.fa", test_dir);
        let file4_path = format!("{}/file4Q.fa", test_dir);
        let file5_path = format!("{}/file5Q.fa", test_dir);
        let file6_path = format!("{}/file6Q.fa", test_dir);
        let file7_path = format!("{}/file7Q.fa", test_dir);
        let file8_path = format!("{}/file8Q.fa", test_dir);

        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file2_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file3_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file4_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file5_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file6_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file7_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file8_path).expect("Failed to write to fof.txt");
        }

        for (file_path, (seq_id, ka_value, sequence)) in [
            (&file1_path, ("seq1", 30, "AAAAAAAAAAAAAAAAAAAAAACACAGATCA")),
            (&file2_path, ("seq2", 12, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT")),
            (&file3_path, ("seq3", 1500, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA")),
            (&file4_path, ("seq4", 1000, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG")),
            (&file5_path, ("seq5", 450, "AAAAAAAAAAAAAAAAAAAAACACCCCTGGG")),
            (&file6_path, ("seq6", 4, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA")),
            (&file7_path, ("seq7", 45110, "AAAAAAAAAAAAAAAAAAAAACACCCCTGGG")),
            (&file8_path, ("seq8", 75, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA")),
        ] {
            let mut file = File::create(file_path).expect("Failed to create FASTA file");
            writeln!(file, ">{} ka:f:{}", seq_id, ka_value).expect("Failed to write header");
            writeln!(file, "{}", sequence).expect("Failed to write sequence");
        }

        let k = 31;
        let m = 15;
        let bf_size = 8;
        let partition_number = 2;
        //let color_number = 6; 
        let color_number = 8;
        let abundance_number = 255;
        let abundance_max = 65535;
        let dense_option = false;
        let threshold = color_number;

        let (_file_paths, index_dir) = build_index(
            vec![
                file1_path.clone(),
                file2_path.clone(),
                file3_path.clone(),
                file4_path.clone(),
                file5_path.clone(),
                file6_path.clone(),
                file7_path.clone(),
                file8_path.clone(),
            ],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
             false,
        )
        .expect("Failed to build index");
        let query_results_path = format!("{}/query_results.csv", index_dir);

        query_index(&file1_path, test_dir, &query_results_path, false, false, 0.5)
            .expect("Failed to query sequences");

        let mut reader =
            csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));

        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let seq_id = record[0].to_string();
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((seq_id, color, abundance));
        }

        assert!(
            !results.is_empty(),
            "Empty results"
        );
        let mut expected_results = vec![
            (">seq1 ka:f:30".to_string(), 0, 29),
        ];
        results.sort();
        expected_results.sort();
        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }
        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");

    }

    /*
    #[test]
    fn test_split_fof_large() -> std::io::Result<()> {
        use std::fs::{self, File};
        use std::io::{BufReader, BufRead, Write};
    
        let test_dir = "test_files_fofs_large";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");
    
        let fof_path = format!("{}/fof.txt", test_dir);
    
        // Generate 1100 file paths
        let mut file_paths: Vec<String> = Vec::new();
        for i in 1..=1100 {
            let file_path = format!("{}/file{}.fa", test_dir, i);
            file_paths.push(file_path);
        }
    
        // Write all file paths to fof.txt
        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            for file_path in &file_paths {
                writeln!(fof_file, "{}", file_path).expect("Failed to write to fof.txt");
            }
        }
    
        let file = File::open(&fof_path)?;
        let reader = BufReader::new(file);
        let lines: Vec<String> = reader.lines().collect::<Result<_, _>>()?;
    
        let result = split_fof(&lines);
    
        assert!(result.is_ok(), "split_fof returned an error");
    
        let (fof_chunks, _) = result.unwrap();
    
        // Create expected chunks
        let expected_chunks = vec![
            file_paths[0..1000].to_vec(),
            file_paths[1000..].to_vec(),
        ];
    
        assert_eq!(
            fof_chunks, expected_chunks,
            "Mismatch in chunk distribution: expected {:?}, got {:?}",
            expected_chunks, fof_chunks
        );
    
        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");
    
        Ok(())
    }
    */


    //#[ignore]
    #[test]
    fn test_build_and_query_index_with_chunks_and_merge2() {
        let test_dir = "test_files_bq_merge2";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);
        let file2_path = format!("{}/file2Q.fa", test_dir);
        let file3_path = format!("{}/file3Q.fa", test_dir);
        let file4_path = format!("{}/file4Q.fa", test_dir);
        let file5_path = format!("{}/file5Q.fa", test_dir);
        let file6_path = format!("{}/file6Q.fa", test_dir);


        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file2_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file3_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file4_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file5_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file6_path).expect("Failed to write to fof.txt");

        }

        for (file_path, (seq_id, ka_value, sequence)) in [
            (&file1_path, ("seq1", 30, "AAAAAAAAAAAAAAAAAAAAAACACAGATCA")),
            (&file2_path, ("seq2", 12, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT")),
            (&file3_path, ("seq3", 1500, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA")),
            (&file4_path, ("seq4", 1000, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG")),
            (&file5_path, ("seq5", 450, "AAAAAAAAAAAAAAAAAAAAACACCCCTGGG")),
            (&file6_path, ("seq6", 4, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA")),

        ] {
            let mut file = File::create(file_path).expect("Failed to create FASTA file");
            writeln!(file, ">{} ka:f:{}", seq_id, ka_value).expect("Failed to write header");
            writeln!(file, "{}", sequence).expect("Failed to write sequence");
        }

        let k = 31;
        let m = 15;
        let bf_size = 8;
        let partition_number = 2;
        let color_number = 6; 
        let abundance_number = 255;
        let abundance_max = 65535;
        let dense_option = false;
        let threshold = color_number;

        let (_file_paths, index_dir) = build_index(
            vec![
                file1_path.clone(),
                file2_path.clone(),
                file3_path.clone(),
                file4_path.clone(),
                file5_path.clone(),
                file6_path.clone(),
       
            ],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");
        let query_results_path = format!("{}/query_results.csv", index_dir);
        query_index(&file6_path, &test_dir, &query_results_path, false, false, 0.5)
            .expect("Failed to query sequences");

       
        let mut reader =
            csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));

        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let seq_id = record[0].to_string();
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((seq_id, color, abundance));
        }

        assert!(
            !results.is_empty(),
            "Empty results"
        );
        let mut expected_results = vec![
            (">seq6 ka:f:4".to_string(), 2, 1476),
            (">seq6 ka:f:4".to_string(), 5, 4),
        ];
        results.sort();
        expected_results.sort();
        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }
        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");

    }


    //#[ignore]
    #[test]
    fn test_build_and_query_index_with_chunks_and_merge3() {
        let test_dir = "test_files_bq_merge3";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);
        let file2_path = format!("{}/file2Q.fa", test_dir);
        let file3_path = format!("{}/file3Q.fa", test_dir);
        let file4_path = format!("{}/file4Q.fa", test_dir);
        let file5_path = format!("{}/file5Q.fa", test_dir);
        let file6_path = format!("{}/file6Q.fa", test_dir);
        let file7_path = format!("{}/file7Q.fa", test_dir);
        let file8_path = format!("{}/file8Q.fa", test_dir);
        let file9_path = format!("{}/file9Q.fa", test_dir);


        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file2_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file3_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file4_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file5_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file6_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file7_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file8_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file9_path).expect("Failed to write to fof.txt");

        }

        for (file_path, (seq_id, ka_value, sequence)) in [
            (&file1_path, ("seq1", 30, "AAAAAAAAAAAAAAAAAAAAAACACAGATCA")),
            (&file2_path, ("seq2", 12, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT")),
            (&file3_path, ("seq3", 1500, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA")),
            (&file4_path, ("seq4", 1000, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG")),
            (&file5_path, ("seq5", 60_000, "TAAAAAAAAAAAAAAAAAAAACACCCCTGGG")),
            (&file6_path, ("seq6", 4, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA")),
            (&file7_path, ("seq7", 1000, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG")),
            (&file8_path, ("seq8", 450, "AAAAAAAAAAAAAAAAAAAAACACCCCTGGG")),
            (&file9_path, ("seq9", 4, "GGGGAAAAAAAAAAAAAAAAAACAAAAAGAA")),

        ] {
            let mut file = File::create(file_path).expect("Failed to create FASTA file");
            writeln!(file, ">{} ka:f:{}", seq_id, ka_value).expect("Failed to write header");
            writeln!(file, "{}", sequence).expect("Failed to write sequence");
        }

        let k = 31;
        let m = 15;
        let bf_size = 4;
        let partition_number = 2;
        let color_number = 9; 
        let abundance_number = 255;
        let abundance_max = 65535;
        let dense_option = false;
        let threshold = color_number;

        let (_file_paths, index_dir) = build_index(
            vec![
                file1_path.clone(),
                file2_path.clone(),
                file3_path.clone(),
                file4_path.clone(),
                file5_path.clone(),
                file6_path.clone(),
                file7_path.clone(),
                file8_path.clone(),
                file9_path.clone(),
       
            ],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");
        let query_results_path = format!("{}/query_results.csv", index_dir);

        query_index(&file5_path, test_dir, &query_results_path, false, false, 0.5)
            .expect("Failed to query sequences");

        let mut reader =
            csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));

        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let seq_id = record[0].to_string();
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((seq_id, color, abundance));
        }

        assert!(
            !results.is_empty(),
            "Empty results"
        );
        let mut expected_results = vec![
            (">seq5 ka:f:60000".to_string(), 4, 59149),
        ];
        results.sort();
        expected_results.sort();
        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }

        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");

    }



    //#[ignore]
    #[test]
    fn test_build_and_query_index_with_chunks_and_merge4() {
        let test_dir = "test_files_bq_merge4";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);
        let file2_path = format!("{}/file2Q.fa", test_dir);
        let file3_path = format!("{}/file3Q.fa", test_dir);
        let file4_path = format!("{}/file4Q.fa", test_dir);
        let file5_path = format!("{}/file5Q.fa", test_dir);



        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file2_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file3_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file4_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file5_path).expect("Failed to write to fof.txt");

        }

        for (file_path, sequences) in [
        (
            &file1_path,
            vec![("seq1", 30, "AAAAAAAAAAAAAAAAAAAAAACACAGATCA")],
        ),
        (
            &file2_path,
            vec![("seq2", 12, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT")],
        ),
        (
            &file3_path,
            vec![("seq3", 1500, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA")],
        ),
        (
            &file4_path,
            vec![("seq4", 1000, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG")],
        ),
        (
            &file5_path,
            vec![
                ("seq5", 60_000, "TAAAAAAAAAAAAAAAAAAAACACCCCTGGG"),
                ("seq6", 4, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA"),
            ],
        ),
        ] {
            let mut file = File::create(file_path).expect("Failed to create FASTA file");
            for (seq_id, ka_value, sequence) in sequences {
                writeln!(file, ">{} ka:f:{}", seq_id, ka_value).expect("Failed to write header");
                writeln!(file, "{}", sequence).expect("Failed to write sequence");
            }
        }


        let k = 31;
        let m = 15;
        let bf_size = 256*256;
        let partition_number = 2;
        let color_number = 5; 
        let abundance_number = 255;
        let abundance_max = 65535;
        let dense_option = false;
        let threshold = color_number;

        let (_file_paths, index_dir) = build_index(
            vec![
                file1_path.clone(),
                file2_path.clone(),
                file3_path.clone(),
                file4_path.clone(),
                file5_path.clone(),
       
            ],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");
        let query_results_path = format!("{}/query_results.csv", index_dir);
        query_index(&file5_path, test_dir, &query_results_path, false, false, 0.5)
            .expect("Failed to query sequences");

        
        let mut reader =
            csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));

        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let seq_id = record[0].to_string();
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((seq_id, color, abundance));
        }

        assert!(
            !results.is_empty(),
            "Empty results"
        );
        let mut expected_results = vec![
            (">seq5 ka:f:60000".to_string(), 4, 59149),
            (">seq6 ka:f:4".to_string(), 2, 1476),
            (">seq6 ka:f:4".to_string(), 4, 4),
        ];
        results.sort();
        expected_results.sort();
        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }


        
        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");

    }




    #[test]
    fn test_merge_partition_bloom_filters() {
        use roaring::RoaringBitmap;
        use std::fs::{create_dir_all, File};
        use std::io::Write;

        let test_dir = "test_merge_partition_bloom_filters";
        create_dir_all(test_dir).expect("Failed to create test directory");

        let chunk1_path = format!("{}/chunk1_p0.txt", test_dir);
        let chunk2_path = format!("{}/chunk2_p0.txt", test_dir);
        let chunk3_path = format!("{}/chunk3_p0.txt", test_dir);


        let partition_idx = 0;
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

        let chunk_files = vec![chunk1_path.clone(), chunk2_path.clone(), chunk3_path.clone()];
        let color_counts = vec![2, 2, 1]; //nb of colors for each chunk

        for (chunk_path, (bf, colors)) in chunk_files.iter().zip(vec![(bf1, 2), (bf2, 2), (bf3, 1)]) {
            let mut file = File::create(chunk_path).expect("Failed to create test chunk file");
            file.write_all(&(colors as u64).to_le_bytes()).expect("Failed to write color count");
            bf.serialize_into(&mut file).expect("Failed to serialize Bloom filter");
        }

     

        let mut merged_bf =RoaringBitmap::new();
    
        // Call the modified function
        merge_partition_bloom_filters(
            chunk_files.clone(),
            partition_idx,
            partitioned_bf_size,
            abundance_number,
            &color_counts,
            &mut merged_bf,
            "test_merge_partition_bloom_filters",
            5
        )
        .expect("Failed to merge partition Bloom filters");
    
        // check the merged Bloom filter
        // let merged_bf = merged_bloom_filter.lock().unwrap();

        //let final_output_path = format!("{}/partition_bloom_filters_p{}.txt", output_dir, partition_idx);
        //let (merged_bf, color_nb) = load_bloom_filter(&final_output_path).expect("Failed to load merged Bloom filter");

        //assert_eq!(
        //    color_nb, 5,
        //    "Expected the merged Bloom filter to represent 5 colors, but got {}",
        //    color_nb
        //);

        let expected_elements: Vec<u32> = vec![1, 2, 9,10,29];
        for elem in expected_elements {
            assert!(
                merged_bf.contains(elem),
                "Expected element {} not found in the merged Bloom filter",
                elem
            );
        }

        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");
    }
   
 
    //#[ignore] 
    #[test]
    fn test_build_and_query_index_with_chunks_and_merge_longseq() {
        let test_dir = "test_files_bq_merge3_ls";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fof_path = format!("{}/fof.txt", test_dir);
        let file1_path = format!("{}/file1Q.fa", test_dir);
        let file2_path = format!("{}/file2Q.fa", test_dir);
        let file3_path = format!("{}/file3Q.fa", test_dir);
        let file4_path = format!("{}/file4Q.fa", test_dir);
        let file5_path = format!("{}/file5Q.fa", test_dir);
        let file6_path = format!("{}/file6Q.fa", test_dir);
        let file7_path = format!("{}/file7Q.fa", test_dir);
        let file8_path = format!("{}/file8Q.fa", test_dir);
        let file9_path = format!("{}/file9Q.fa", test_dir);


        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            writeln!(fof_file, "{}", file1_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file2_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file3_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file4_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file5_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file6_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file7_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file8_path).expect("Failed to write to fof.txt");
            writeln!(fof_file, "{}", file9_path).expect("Failed to write to fof.txt");

        }

        for (file_path, (seq_id, ka_value, sequence)) in [
            (&file1_path, ("seq1", 30,    "AAAAAAAAAAAAAAAAAAAAAACACAGATTT")),
            (&file2_path, ("seq2", 12,    "AAAAAAAAAAAAAAAAAAAAACACAGATCAT")),
            (&file3_path, ("seq3", 1500,  "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA")),
            (&file4_path, ("seq4", 1000,  "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG")),
            (&file5_path, ("seq5", 60_000,"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")),
            (&file6_path, ("seq6", 4,     "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA")),
            (&file7_path, ("seq7", 1000,  "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG")),
            (&file8_path, ("seq8", 450,   "AAAAAAAAAAAAAAAAAAAAACACCCCTGGG")),
            (&file9_path, ("seq9", 4,     "GGGGAAAAAAAAAAAAAAAAAACAAAAAGAA")),

        ] {
            let mut file = File::create(file_path).expect("Failed to create FASTA file");
            writeln!(file, ">{} ka:f:{}", seq_id, ka_value).expect("Failed to write header");
            writeln!(file, "{}", sequence).expect("Failed to write sequence");
        }

        let k = 31;
        let m = 15;
        let bf_size = 32;
        let partition_number = 2;
        let color_number = 9; 
        let abundance_number = 255;
        let abundance_max = 65535;
        let dense_option = false;
        let threshold = color_number;

        let (_file_paths, index_dir) = build_index(
            vec![
                file1_path.clone(),
                file2_path.clone(),
                file3_path.clone(),
                file4_path.clone(),
                file5_path.clone(),
                file6_path.clone(),
                file7_path.clone(),
                file8_path.clone(),
                file9_path.clone(),
       
            ],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");
        let query_results_path = format!("{}/query_results.csv", index_dir);

        query_index(&file2_path, &test_dir, &query_results_path, false, false, 0.5)
            .expect("Failed to query sequences");

        let mut reader =
            csv::Reader::from_reader(File::open(&query_results_path).expect("Failed to open query results"));

        let mut results: Vec<(String, usize, usize)> = Vec::new();
        for record in reader.records() {
            let record = record.expect("Failed to read record");
            let seq_id = record[0].to_string();
            let color: usize = record[1].parse().expect("Failed to parse color");
            let abundance: usize = record[2].parse().expect("Failed to parse abundance");
            results.push((seq_id, color, abundance));
        }

        assert!(
            !results.is_empty(),
            "Empty results"
        );
        let expected_results = vec![
            //("seq5".to_string(), 4, 57549),
            (">seq2 ka:f:12".to_string(), 1, 12),
        ];

        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }


        
        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");

    }


    #[test]
    fn test_color_graph() {
        let test_dir = "test_color_graph";
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fasta_path = format!("{}/test.fasta", test_dir);
        let output_path = format!("{}/colored_graph_output.fasta", test_dir);

        {
            let mut file = File::create(&fasta_path).expect("Failed to create test.fasta");
            writeln!(file, ">seq1 ka:f:29").expect("Failed to write header");
            writeln!(file, "TAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACAGATCA").expect("Failed to write sequence");
            writeln!(file, ">seq2 ka:f:250").expect("Failed to write header");
            writeln!(file, "TTTTTAATGATCGATTTTTTTTTTTACCCCTGG").expect("Failed to write sequence");
        }

        let k = 31;
        let m = 15;
        let bf_size = 1024 * 1024;
        let partition_number = 4;
        let color_number = 1;
        let abundance_number = 255;
        let abundance_max = 65535;
        let _batch_size = 2;
        let dense_option = false;
        let threshold = color_number;


        let (_file_paths, _index_dir) = build_index(
            vec![
                fasta_path.clone(),
            ],
            k,
            m,
            bf_size,
            partition_number,
            color_number,
            abundance_number,
            abundance_max,
            test_dir,
            dense_option,
            threshold,
            false,
        )
        .expect("Failed to build index");
        query_index(
            &fasta_path,
            test_dir,
            &output_path,
            true, //  graph coloring
            false, 
            0.5,
        )
        .expect("Failed to color graph");

        let expected_output = vec![
           ">seq1 ka:f:29 col:0:29", "TAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACAGATCA", ">seq2 ka:f:250 col:0:249", "TTTTTAATGATCGATTTTTTTTTTTACCCCTGG",
        ];

        let mut actual_output = Vec::new();
        let output_file = File::open(&output_path).expect("Failed to open output file");
        let reader = BufReader::new(output_file);

        for line in reader.lines() {
            actual_output.push(line.expect("Failed to read line"));
        }

        assert_eq!(
            expected_output, actual_output,
            "Mismatch between expected and actual output"
        );

        fs::remove_dir_all(test_dir).expect("Failed to clean up test directory");
    }
    

    // use csv::Reader;
    // #[test]
    // fn test_merge_multiple_indexes_from_fof() -> io::Result<()> {
    //     let base_dir = "test_merge_indexes";
    //     let index1_files_dir = format!("{}/index1_files", base_dir);
    //     let index2_files_dir = format!("{}/index2_files", base_dir);
    //     let index1_index_dir = format!("{}/index1_index", base_dir);
    //     let index2_index_dir = format!("{}/index2_index", base_dir);
    //     let merged_index_dir = format!("{}/merged_index", base_dir);
    //     let indexes_fof = format!("{}/indexes.txt", base_dir);

    //     fs::create_dir_all(&index1_files_dir)?;
    //     fs::create_dir_all(&index2_files_dir)?;
    //     fs::create_dir_all(&index1_index_dir)?;
    //     fs::create_dir_all(&index2_index_dir)?;

    //     let file1_path = format!("{}/file1.fa", index1_files_dir);
    //     {
    //         let mut file = File::create(&file1_path)?;
    //         writeln!(file, ">seq1 ka:f:30")?;
    //         writeln!(file, "ACGTACG")?;
    //     }

    //     let file2_path = format!("{}/file2.fa", index2_files_dir);
    //     let file3_path = format!("{}/file3.fa", index2_files_dir);
    //     {
    //         let mut file = File::create(&file2_path)?;
    //         writeln!(file, ">seq2 ka:f:30")?;
    //         writeln!(file, "ACGTACG")?;
    //     }
    //     {
    //         let mut file = File::create(&file3_path)?;
    //         writeln!(file, ">seq3 ka:f:30")?;
    //         writeln!(file, "ACGTACG")?;
    //     }

    //     let k = 7;
    //     let m = 3;
    //     let bf_size = 1024;         
    //     let partitions = 2;         
    //     let abundance = 255;
    //     let abundance_max = 65535;
    //     let dense_option = false;

    //     // for index1, color count = 1
    //     let index1_file_paths = vec![file1_path.clone()];
    //     let (_dummy, index_dir1) = build_index(
    //         index1_file_paths,
    //         k,
    //         m,
    //         bf_size,
    //         partitions,
    //         1, // color_number = 1
    //         abundance,
    //         abundance_max,
    //         &String::from(index1_index_dir),
    //         dense_option,
    //         1
    //     )?;
    //     //index2
    //     let index2_file_paths = vec![file2_path.clone(), file3_path.clone()];
    //     let (_dummy, index_dir2) = build_index(
    //         index2_file_paths,
    //         k,
    //         m,
    //         bf_size,
    //         partitions,
    //         2, // color_number = 2
    //         abundance,
    //         abundance_max,
    //         &String::from(index2_index_dir),
    //         dense_option,
    //         2
    //     )?;

    //     {
    //         let mut fof = File::create(&indexes_fof)?;
    //         writeln!(fof, "{}", index_dir1)?;
    //         writeln!(fof, "{}", index_dir2)?;
    //     }

    //     merge_all_partitions(&indexes_fof, Some(&merged_index_dir))?;

    //     let merged_csv_path: String = format!("{}/index_info.csv", merged_index_dir);
    //     let mut rdr = Reader::from_reader(File::open(&merged_csv_path)?);
    //     let record = rdr
    //         .records()
    //         .next()
    //         .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Merged index metadata CSV is empty"))??;
    //     // in order: k, m, bf_size, partition_number, color_number, abundance_number, abundance_max
    //     let merged_color: usize = record.get(4).unwrap().parse().unwrap();
    //     assert_eq!(merged_color, 3, "Expected merged color count to be 3, got {}", merged_color);

    //     // check that each partition file in the merged index exists
    //     for partition in 0..partitions {
    //         let part_path = format!("{}/partition_bloom_filters_p{}.txt", merged_index_dir, partition);
    //         assert!(Path::new(&part_path).exists(), "Merged partition file {} does not exist", part_path);
    //     }

    //     fs::remove_dir_all(base_dir)?;

    //     Ok(())
    // }

}
