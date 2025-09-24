use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::Path;

use bio::io::fasta;
use rayon::prelude::*;

use crate::reindeer2::{
    approximate_value, compute_base, compute_base_position, kmer_minimizers_seq_level,
    load_bloom_filter, load_dense_index, process_fasta_in_batches, read_file,
    update_color_abundances,
};
use crate::OutputFormat;

// === QUERY ===

// --- MAIN FUNCTION ---

// --- QUERY FUNCTIONS ---

/// Builds a map from partition_index -> Vec of (sequence_id, kmer_hash).
fn build_partitions_kmers(
    batch: Vec<fasta::Record>,
    k: usize,
    m: usize,
    partition_number: u64,
) -> HashMap<usize, Vec<(String, u64)>> {
    let mut partition_kmers: HashMap<usize, Vec<(String, u64)>> = HashMap::new();

    // Build all partition-kmers upfront
    for record in batch {
        let id = record.id();
        let desc = record.desc().unwrap_or("");
        // Build the header string only once per record
        let full_header = format!(">{} {}", id, desc).trim().to_string();
        // let full_header = format!(">{}", id).trim().to_string();

        let seq_str = std::str::from_utf8(record.seq())
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid UTF-8 sequence"))
            .unwrap(); // or handle error gracefully

        // k-mer hash + minimizer
        //let nt_hash_iterator = NtHashIterator::new(seq_str.as_bytes(), k).unwrap();
        //let min_iter = MinimizerBuilder::<u64>::new()
        //    .minimizer_size(m)
        //    .width((k - m + 1).try_into().unwrap())
        //    .iter(seq_str.as_bytes());

        for (kmer_hash, minimizer) in kmer_minimizers_seq_level(seq_str.as_bytes(), k, m) {
            //for (kmer_hash, (minimizer, _position)) in nt_hash_iterator.zip(min_iter) {
            let partition_index = (minimizer % partition_number) as usize;
            partition_kmers
                .entry(partition_index)
                .or_default()
                .push((full_header.clone(), kmer_hash));
        }
    }
    partition_kmers
}

// TODO rename
fn fold_into_hashmap(
    mut local_results: HashMap<String, Vec<Vec<u16>>>,
    partition_index: usize,
    kmers: Vec<(String, u64)>,
    base: f64,
    bf_dir: &str,
    bf_size: u64,
    partition_number: usize,
    color_number: usize,
    abundance_number: usize,
    is_dense: bool,
) -> HashMap<String, Vec<Vec<u16>>> {
    // Load the partition's Bloom filter
    let path_bf = format!(
        "{}/partition_bloom_filters_p{}.txt",
        bf_dir, partition_index
    );
    let maybe_bf = load_bloom_filter(&path_bf);

    if let Ok((bitmap, _maybe_aux_data)) = maybe_bf {
        let hashmap: HashMap<u64, Box<Vec<u8>>> = if is_dense {
            let path_dense_index =
                format!("{}/partition_dense_index_p{}.txt", bf_dir, partition_index);
            load_dense_index(&path_dense_index).expect(&format!(
                "Failed to load dense index for partition {}",
                partition_index
            ))
        } else {
            HashMap::new()
        };
        //  For each k-mer in this partition
        for (sequence_id, kmer_hash) in kmers {
            let color_abundances = if hashmap.contains_key(&kmer_hash) {
                let log_abundance_vector =
                    hashmap.get(&kmer_hash).expect("failed to read the hashmap");
                let mut color_abundances = vec![Vec::new(); color_number];
                for (color, log_abundance) in log_abundance_vector.iter().enumerate() {
                    if *log_abundance == 0 {
                        // TODO discuss weird value
                        color_abundances[color].push(666);
                    } else {
                        color_abundances[color].push((*log_abundance - 1) as usize);
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
                        .map(|log_abund| {
                            if log_abund == 666 {
                                0
                            } else {
                                approximate_value(log_abund, base)
                            }
                        })
                        .collect()
                })
                .collect();

            // Accumulate results in local_results
            let entry = local_results
                .entry(sequence_id.clone())
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
    let mut file = File::open(Path::new(dir_path).join("kmer_counts_per_color.txt"))?;

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

/// Write abundances per kmer like RD1.
fn write_raw_abundance(
    bf_dir: &str,
    sequence_results: &HashMap<String, Vec<Vec<u16>>>,
    writer: &mut impl Write,
) -> io::Result<()> {
    // TODO make RD1 and normalize mutually exclusive at compile time
    // we need the count of kmers for an outpout like the one of Reindeer 1
    load_kmer_counts_vector(bf_dir).expect("Failed to load from disk the kmer counts vector");

    for (seq_header, color_vectors) in sequence_results {
        write!(writer, "{seq_header}")?;
        for abund_values in color_vectors.iter() {
            // new color => a space
            write!(writer, " ")?;
            let mut start = 0;
            let mut current = abund_values[0];

            for i in 1..=abund_values.len() {
                if i == abund_values.len() || abund_values[i] != current {
                    let val_str = if current == 0 {
                        "*"
                    } else {
                        &current.to_string()
                    };
                    if start + 1 == i {
                        write!(writer, "{}:{}", start, val_str)?;
                    } else {
                        write!(writer, "{}-{}:{}", start, i - 1, val_str)?;
                    }
                    if i < abund_values.len() {
                        write!(writer, ",")?;
                        start = i;
                        current = abund_values[i];
                    }
                }
            }
        }
        writeln!(writer)?;
    }
    Ok(())
}

/// Write abundances per kmer like RD1.
fn write_median_abundance(
    bf_dir: &str,
    sequence_results: &HashMap<String, Vec<Vec<u16>>>,
    color_number: usize,
    normalize: bool,
    coverage: f32,
    writer: &mut impl Write,
) -> io::Result<()> {
    writeln!(writer, "header,file,abundance")?;
    // we need the count of kmers if we want to normalize them
    let kmer_counts = if normalize {
        load_kmer_counts_vector(bf_dir).expect("Failed to load from disk the kmer counts vector")
    } else {
        vec![color_number, 0] // TODO bizarre
    };
    // Compute medians for each sequence and each color, then write them out
    for (seq_header, color_vectors) in sequence_results {
        for (color_idx, abund_values) in color_vectors.iter().enumerate() {
            if !abund_values.is_empty() {
                let mut zeros_count = 0;
                let mut non_zero_values: Vec<u16> = Vec::new();
                abund_values.iter().for_each(|value| {
                    if *value == 0 {
                        zeros_count += 1;
                    } else {
                        non_zero_values.push(*value);
                    }
                });
                if !non_zero_values.is_empty()
                    && (((zeros_count as f32) / (abund_values.len() as f32)) < coverage)
                {
                    let mut abund_sorted = non_zero_values.clone();
                    abund_sorted.sort_unstable();
                    let median = if abund_sorted.iter().all(|&x| x == 0) {
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
                        let median = if normalize {
                            median as f64 / kmer_counts[color_idx] as f64 * 1_000_000 as f64
                        } else {
                            median as f64
                        };
                        writeln!(writer, "{},{},{}", seq_header, color_idx, median)?;
                    }
                }
            }
        }
    }
    Ok(())
}

fn write_kmer_query(
    fasta_file: &str,
    bf_dir: &str,
    color_number: usize,
    batch_size: usize,
    output_format: OutputFormat,
    coverage: f32,
    sequence_results: HashMap<String, Vec<Vec<u16>>>,
    mut writer: &mut impl Write,
) -> io::Result<()> {
    match output_format {
        OutputFormat::Colored => {
            // TODO lrobidou discuss what this comment means:
            // If your graph coloring wants to read from `sequence_results`:
            // Flush the writer to separate batch outputs if needed
            graph_coloring(fasta_file, batch_size, writer, &sequence_results)
        }
        OutputFormat::Raw => write_raw_abundance(bf_dir, &sequence_results, &mut writer),
        OutputFormat::Median | OutputFormat::NormalizedMedian => {
            let normalize = output_format == OutputFormat::NormalizedMedian;
            write_median_abundance(
                bf_dir,
                &sequence_results,
                color_number,
                normalize,
                coverage,
                &mut writer,
            )
        }
    }?;
    writer.flush()?;
    Ok(())
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
) -> io::Result<()> {
    let reader = read_file(fasta_file)?;
    let mut writer = BufWriter::new(File::create(output_file)?);
    // write the headr of the output csv file
    let base = compute_base(abundance_number, abundance_max);

    // Process FASTA in chunks of `batch_size` records
    process_fasta_in_batches(reader, batch_size, |batch| {
        let partition_kmers = build_partitions_kmers(batch, k, m, partition_number as u64);

        // --- PARALLEL PHASE: process each partition's k-mers in parallel ---
        // This will store final results for *all* sequences in this batch.
        // Key: sequence header; Value: vector of vector of abundances
        let sequence_results = partition_kmers
            .into_par_iter()
            // 1) Create a local HashMap in each thread
            .fold(
                HashMap::<String, Vec<Vec<u16>>>::new,
                |local_results, (partition_index, kmers)| {
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
            .reduce(HashMap::<String, Vec<Vec<u16>>>::new, merge_results);

        // Now `sequence_results` has the combined data for this batch.
        // Let's compute the output in the requested format.
        write_kmer_query(
            fasta_file,
            bf_dir,
            color_number,
            batch_size,
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
    mut acc: HashMap<String, Vec<Vec<u16>>>,
    local: HashMap<String, Vec<Vec<u16>>>,
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
    writer: &mut impl Write,
    sequence_results: &HashMap<String, Vec<Vec<u16>>>,
) -> std::io::Result<()> {
    let reader = read_file(fasta_file)?; // your existing read_file

    process_fasta_in_batches(reader, batch_size, |batch| {
        for record in batch {
            let id = record.id();
            let desc = record.desc().unwrap_or("");
            let full_header = format!(">{} {}", id, desc).trim().to_string();
            let seq_str = std::str::from_utf8(record.seq()).expect("Invalid UTF-8 sequence");

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
                    let median = if abund_sorted.iter().all(|&x| x == 0) {
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
    Ok(())
}
