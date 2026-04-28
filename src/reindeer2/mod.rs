mod dense_index;
mod filter;
mod index;
mod merge;
mod minimizer_iter;
mod query;
mod test_utils;
use bio::io::fasta;
use flate2::read::GzDecoder;
use itertools::Itertools;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::num::NonZero;
use std::panic;
use std::path::Path;
use std::sync::{atomic, Arc, Mutex};

pub use merge::merge_multiple_indexes;

#[cfg(any(debug_assertions, test))]
use thousands::Separable;

use zstd::stream::decode_all;

use crate::reindeer2::dense_index::DenseIndex;
use crate::reindeer2::filter::Filters;
use crate::reindeer2::minimizer_iter::{KmerSampler, MinimizerSampler, NoSampler, Sampler};
use crate::reindeer2::query::{fimpera, load_kmer_counts_vector, ApproxAbundance};

const NB_FILE_IN_AN_INDEX: usize = 1024;

// TODO move to another place ?
fn create_and_reserve_tar_get_file<P>(path: P, nb_object: u64) -> File
where
    P: AsRef<Path>,
{
    let file = OpenOptions::new()
        .read(true)
        .write(true)
        .create_new(true)
        .open(&path)
        .unwrap_or_else(|_| panic!("should have been able to create file {}, maybe it already exists or you don't have permissions to create it", path.as_ref().display()));
    let mut writer = BufWriter::new(file);
    tar_get::reserve_capacity(&mut writer, nb_object)
        .expect("should have been able to write on disk");

    writer.into_inner().unwrap_or_else(|_| {
        panic!(
            "should have been able to flush the iternal buffer to {}",
            path.as_ref().display()
        )
    })
}

#[derive(Clone, Debug, PartialEq)]
/// Breakpoints reports position in which the abundance varies along the queried k-mers.
/// Normalize prints the values of k-mers.
///
/// They are not compatible in the current implementation.
pub enum BreakpointsNormalize {
    /// Ask to report the position at which the abundance varies along the queried k-mers.
    /// Parameter: the parameter given to the PELT algorithm to detect breakpoints.
    Breakpoints(f64),
    /// Normalize the output by a factor decided by the parameter and the number of kindexed k-mers.
    Normalize(u64),
}

/// Matrix formats supported by REINDEER 2.
#[derive(Clone, Debug, PartialEq)]
pub enum MatrixFormat {
    /// Raw format: the abundance of each k-mer is present in the output.
    Raw(Option<BreakpointsNormalize>),
    /// Average format: the average of the k-mers of the sequence is outputted.
    Average {
        /// This format can be normalized, i.e. the abundance is adjusted
        /// by a factor depending on a constant and on the number of k-mers in the indexed dataset.
        normalized: Option<u64>,
    },
    /// Median format: the median of the k-mers of the sequence is outputted.
    Median {
        /// This format can be normalized, i.e. the abundance is adjusted
        /// by a factor depending on a constant and on the number of k-mers in the indexed dataset.
        normalized: Option<u64>,
    },
}

/// Output formats supported by REINDEER 2.
#[derive(Clone, Debug, PartialEq)]
pub enum OutputFormat {
    /// Colored format: annotate the input file with abundances rather than producing the standard output file.
    Colored {
        /// This format can be normalized, i.e. the abundance is adjusted
        /// by a factor depending on a constant and on the number of k-mers in the indexed dataset.
        normalized: Option<u64>,
    },
    /// Median format: for each color, returns the median of k-mer abundance per read.
    Median {
        /// This format can be normalized, i.e. the abundance is adjusted
        /// by a factor depending on a constant and on the number of k-mers in the indexed dataset.
        normalized: Option<u64>,
    },
    /// Matrix formats: output a matrix of k-mers. Row: queried sequences, columns: indexed datasets.
    AbundanceMatrix {
        /// The format of the matrix.
        format: MatrixFormat,
    },
}

/// REINDEER 2's  parameters
#[derive(Serialize, Deserialize, Clone, Debug)]
#[cfg_attr(any(test, debug_assertions), derive(PartialEq))]
pub struct Parameters {
    /// Size of the Bloom filter in bits
    pub bf_size: u64,
    /// Number of partition in the index
    pub partition_number: usize,
    /// Size of the k-mers
    pub k: usize,
    /// Size of the minimizers
    pub m: usize,
    /// Number of color in the index
    pub nb_color: usize,
    /// Number of abundance levels in the index
    pub abundance_number: NonZero<usize>,
    /// Minimum abundance for a k-mer to be indexed
    pub abundance_min: u16,
    /// Maximal abundance of k-mers. If a k-mer has a greater abundance, its abundance is capped.
    pub abundance_max: NonZero<u16>,
    /// True if and only if the index a dense index
    pub dense_option: bool,
    /// True if and only if the index is buit in canonical mode, i.e. we consider a strand equal to its reverse complement
    pub canonical: bool,
    /// The strategy to sample k-mers or minimizers.
    /// Allows to save query time, indexation time and memory, at the cost of some approxitmation.
    /// If None, no sampling is performed.
    pub sampling_strategy: Option<SamplingStrategy>,
    /// The z parameter of findere.
    /// The index will contains (k-z)-mers, and reconstruct k-mers on the fly.
    /// This allows to decrease the false positive rate of queries.
    pub findere_z: usize,
    pub capacity: usize,
}

/// The strategy to sample k-mers or minimizers.
/// Allows to save query time, indexation time and memory, at the cost of some approxitmation.
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq)]
pub enum SamplingStrategy {
    /// Only take into account minimizers which hash ends by at least `last_bits_to_zero` 0s.
    MinimizerSampling {
        /// The minimum number of 0 a minimizer should have at the end of its hash to be sampled
        last_bits_to_zero: u64,
    },
    /// Only take into account k-mers which hash ends by at least `last_bits_to_zero` 0s.
    KmerSampling {
        /// The minimum number of 0 a k-mer should have at the end of its hash to be sampled
        last_bits_to_zero: u64,
    },
}

impl Parameters {
    /// Constructs a new `Parameters` from its inner values.
    pub const fn new(
        bf_size: u64,
        partition_number: usize,
        k: usize,
        m: usize,
        nb_color: usize,
        abundance_number: NonZero<usize>,
        abundance_min: u16,
        abundance_max: NonZero<u16>,
        // TODO rename is_dense
        dense_option: bool,
        canonical: bool,
        sampling_strategy: Option<SamplingStrategy>,
        findere_z: usize,
        capacity: usize,
    ) -> Self {
        Self {
            bf_size,
            partition_number,
            k,
            m,
            nb_color,
            abundance_number,
            abundance_min,
            abundance_max,
            dense_option,
            canonical,
            sampling_strategy,
            findere_z,
            capacity,
        }
    }

    fn can_merge(&self, parameters2: &Parameters) -> bool {
        self.k == parameters2.k
            && self.m == parameters2.m
            && self.bf_size == parameters2.bf_size
            && self.partition_number == parameters2.partition_number
            && self.abundance_number == parameters2.abundance_number
            && self.abundance_max == parameters2.abundance_max
            && self.dense_option == parameters2.dense_option
            && self.canonical == parameters2.canonical
            && self.sampling_strategy == parameters2.sampling_strategy
    }
}

/// A REINDEER 2 index handle. Allows to build and query the index.
#[derive(Serialize, Deserialize)]
#[cfg_attr(any(test, debug_assertions), derive(PartialEq, Debug))]
pub struct Reindeer2 {
    /// Parameters used to build the index
    parameters: Parameters,
    /// Names of indexed files
    indexed_file_names: Vec<String>,
    /// Directory containing the index
    #[serde(skip)]
    index_dir: String, // TODO use a path ?
}

/// Declares a variable. The variable is declared as mut iif `#[cfg(any(debug_assertions, test))]`.
macro_rules! mut_if_debug {
    ($name:ident = $val:expr) => {
        #[cfg(any(debug_assertions, test))]
        let mut $name = $val;

        #[cfg(not(any(debug_assertions, test)))]
        let $name = $val;
    };
}

/// Infos about a REINDEER 2 index.
///
/// Implements serde::Serialize and serde::Deserialize so it can be e.g. shown to users.
#[derive(Serialize, Deserialize)]
pub struct Infos {
    /// Parameters used to build the index.
    pub parameters: Parameters,
    /// Vector of infos about a file. Each element is a tuple (file name, number of k-mers in the file).
    pub indexed_file_names: Vec<(String, usize)>,
    /// Directory contaning the index.
    pub index_dir: String,
}

use std::fmt;

impl fmt::Display for Infos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "index directory: \"{}\"", self.index_dir)?;
        writeln!(f, "parameters: {:#?}", self.parameters)?;
        writeln!(f, "indexed filenames and k-mers: [")?;
        for (name, size) in &self.indexed_file_names {
            writeln!(f, "  (\"{}\", {}),", name, size)?;
        }
        writeln!(f, "]")
    }
}

impl Reindeer2 {
    /// Constructs a new index handle.
    pub const fn new(parameters: Parameters, index_dir: String) -> Self {
        Self {
            parameters,
            indexed_file_names: Vec::new(),
            index_dir,
        }
    }

    /// Returns the infos of an index.
    pub fn get_index_infos(&self) -> Infos {
        let Self {
            parameters,
            indexed_file_names,
            index_dir,
        } = self;
        let parameters = parameters.clone();
        let index_dir = index_dir.clone();
        let kmer_count = load_kmer_counts_vector(&self.index_dir)
            .expect("Failed to load from disk the kmer counts vector");
        let indexed_file_names = indexed_file_names.iter().cloned().zip(kmer_count).collect();
        Infos {
            parameters,
            indexed_file_names,
            index_dir,
        }
    }

    #[cfg(any(debug_assertions, test))]
    /// Sets the name of indexed file in the index.
    ///
    /// Only available in test and debug mode.
    ///
    /// # Warning
    /// This breaks your ability to reason about the index. Do not use except for unit tests.
    pub fn set_indexed_file_names(&mut self, indexed_file_names: Vec<String>) {
        self.indexed_file_names = indexed_file_names;
    }

    fn index_infos_file(bf_dir: &str) -> String {
        format!("{}/index_info.json", bf_dir)
    }

    /// Loads an index from the JSON file in the directory of the index
    pub fn load_from_disk(bf_dir: impl Into<String>) -> io::Result<Self> {
        let bf_dir = bf_dir.into();
        log::info!("Loading index from disk ({})...", bf_dir);
        let input_path = Self::index_infos_file(&bf_dir);
        let file = File::open(input_path)?;
        let reader = BufReader::new(file);
        let mut rd2: Reindeer2 = serde_json::from_reader(reader)?;
        debug_assert_eq!(rd2.index_dir, "");
        rd2.index_dir = bf_dir;
        log::info!("Index loaded from disk.");
        Ok(rd2)
    }

    /// Save an index metadata into a JSON file in the directory of the index
    fn save_to_disk(&self) -> io::Result<()> {
        log::info!("Saving index information to {}", self.index_dir);
        let output_path = Self::index_infos_file(&self.index_dir);
        let file = File::create(&output_path)?;
        let writer = BufWriter::new(file);
        serde_json::to_writer(writer, self)?;
        log::info!("Index information written to {}", output_path);
        Ok(())
    }

    /// Build the index by indexing the files in `file_paths` using the parameters passed during the construction of the struct.
    pub fn build(
        &mut self,
        file_paths: Vec<String>,
        chunks_size: usize,
        threshold: usize,
        allow_count_right_after_angle_bracket: bool,
    ) -> io::Result<(Vec<String>, String)> {
        mut_if_debug!(total_kmers = atomic::AtomicU64::new(0));
        mut_if_debug!(atomic_dense_kmers_count = atomic::AtomicU64::new(0));
        mut_if_debug!(atomic_sparse_kmers_count = atomic::AtomicU64::new(0));

        let parameters = &self.parameters;

        #[cfg(any(debug_assertions, test))]
        let mut atomic_sparse_one_seen = atomic::AtomicU64::new(0);

        #[cfg(any(debug_assertions, test))]
        let mut atomic_sparse_fp_seen = atomic::AtomicU64::new(0);

        let kmer_counts_vector: Arc<Mutex<Vec<usize>>> =
            Arc::new(Mutex::new(vec![0; parameters.nb_color]));
        let chunks = split_fof(&file_paths, chunks_size)?;
        // the number of color in each chunck
        let color_chunks = chunks.iter().map(Vec::len).collect_vec();
        let base = compute_base(parameters.abundance_number, parameters.abundance_max);

        #[cfg(any(debug_assertions, test))]
        log::debug!("Using log base {}", base);

        let partitioned_bf_size = (parameters.bf_size as usize) / parameters.partition_number;

        #[cfg(any(debug_assertions, test))]
        log::debug!("In debug mode... the tool may take (much) longer than usual.");

        log::info!("Initializing Bloom filter slices...");

        let (_, dir_path) = create_dir_and_files(parameters.partition_number, &self.index_dir)?;

        // Shared data structures protected by Mutex for safe parallel access
        let maybe_dense_indexes: Option<Arc<DenseIndex>> = if parameters.dense_option {
            Some(Arc::new(DenseIndex::with_partition_number(
                parameters.partition_number,
            )))
        } else {
            None
        };

        for (chunk_i, chunk) in chunks.iter().enumerate() {
            // OPTIMIZE we are losing the underlying allocation of the Filters here
            // but currently each chunk could have a different size, so we have to recreate them
            // if this is a bottleneck, we should find a way to reuse the Filters
            let bloom_filters = Filters::with_number_partition(
                parameters.partition_number,
                chunk.len(),
                // TODO unit: is it in bits ?
                // TODO can I use usize here ?
                parameters.bf_size as usize,
                parameters.abundance_number.get(),
            );

            // For each file in this chunk, process in *parallel* (soon)
            // TODO build the appropriate iterator to parallelize (or not) if dense is set
            chunk.par_iter().enumerate().for_each(|(path_num, path)| {
                match std::fs::metadata(path) {
                    Ok(metadata) => {
                        if metadata.is_file() {
                            let reader = match read_file(path) {
                                Ok(r) => r,
                                Err(e) => {
                                    log::error!("Failed to open file {}: {}", path, e);
                                    return;
                                }
                            };
                            let fasta_reader = fasta::Reader::new(reader);
                            let first_record = fasta_reader.records().next();

                            if let Some(Ok(record)) = first_record {
                                let h_type = get_first_header_type(
                                    &record,
                                    allow_count_right_after_angle_bracket,
                                    path,
                                );
                                let (h_type, count_right_after_angle_bracket) = match h_type {
                                    Some(h_type) => h_type,
                                    None => return,
                                };

                                if let Err(e) = match &parameters.sampling_strategy {
                                    None => {
                                        let sampler = NoSampler::new();
                                        index::process_fasta_file::<NoSampler>(
                                            path,
                                            &maybe_dense_indexes,
                                            &bloom_filters,
                                            parameters.k,
                                            parameters.findere_z,
                                            parameters.m,
                                            parameters.partition_number,
                                            parameters.nb_color,
                                            threshold,
                                            parameters.abundance_min,
                                            parameters.abundance_max,
                                            path_num,
                                            path_num + chunk_i * color_chunks[0],
                                            base,
                                            chunk_i,
                                            h_type,
                                            1_000_000, // max size for flushing k-mers to bloom filter
                                            &total_kmers,
                                            &atomic_dense_kmers_count,
                                            &atomic_sparse_kmers_count,
                                            &kmer_counts_vector,
                                            parameters.canonical,
                                            count_right_after_angle_bracket,
                                            &sampler,
                                        )
                                    }
                                    Some(SamplingStrategy::MinimizerSampling {
                                        last_bits_to_zero,
                                    }) => {
                                        let sampler = MinimizerSampler::new(*last_bits_to_zero);
                                        index::process_fasta_file::<MinimizerSampler>(
                                            path,
                                            &maybe_dense_indexes,
                                            &bloom_filters,
                                            parameters.k,
                                            parameters.findere_z,
                                            parameters.m,
                                            parameters.partition_number,
                                            parameters.nb_color,
                                            threshold,
                                            parameters.abundance_min,
                                            parameters.abundance_max,
                                            path_num,
                                            path_num + chunk_i * color_chunks[0],
                                            base,
                                            chunk_i,
                                            h_type,
                                            1_000_000, // max size for flushing k-mers to bloom filter
                                            &total_kmers,
                                            &atomic_dense_kmers_count,
                                            &atomic_sparse_kmers_count,
                                            &kmer_counts_vector,
                                            parameters.canonical,
                                            count_right_after_angle_bracket,
                                            &sampler,
                                        )
                                    }
                                    Some(SamplingStrategy::KmerSampling { last_bits_to_zero }) => {
                                        let sampler = KmerSampler::new(*last_bits_to_zero);
                                        index::process_fasta_file::<KmerSampler>(
                                            path,
                                            &maybe_dense_indexes,
                                            &bloom_filters,
                                            parameters.k,
                                            parameters.findere_z,
                                            parameters.m,
                                            parameters.partition_number,
                                            parameters.nb_color,
                                            threshold,
                                            parameters.abundance_min,
                                            parameters.abundance_max,
                                            path_num,
                                            path_num + chunk_i * color_chunks[0],
                                            base,
                                            chunk_i,
                                            h_type,
                                            1_000_000, // max size for flushing k-mers to bloom filter
                                            &total_kmers,
                                            &atomic_dense_kmers_count,
                                            &atomic_sparse_kmers_count,
                                            &kmer_counts_vector,
                                            parameters.canonical,
                                            count_right_after_angle_bracket,
                                            &sampler,
                                        )
                                    }
                                } {
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
            #[cfg(any(debug_assertions, test))]
            {
                bloom_filters.update_sparse_counts(
                    &atomic_sparse_one_seen,
                    &atomic_sparse_fp_seen,
                    parameters.abundance_number.get(),
                )?;
            }

            if chunks.len() > 1 {
                if let Err(e) = bloom_filters.write_to_disk_as_multiple_small_files(
                    &dir_path,
                    &color_chunks,
                    chunk_i,
                ) {
                    eprintln!("Error writing Bloom filters for chunk {}: {}", chunk_i, e);
                }
            } else {
                // If there is only one chunk, write the final file directly
                if let Err(e) = bloom_filters.write_to_disk_tar_get_files(&dir_path) {
                    eprintln!("Error writing Bloom filters for chunk {}: {}", chunk_i, e);
                }
            }
            log::trace!("Chunk {} done", chunk_i);
        }

        // After processing all chunks, write the dense indexes to disk
        if let Some(dense_indexes) = maybe_dense_indexes {
            dense_indexes.write_to_disk(&dir_path)?;
        }

        write_kmer_counts_to_disk(&dir_path, &kmer_counts_vector)?;

        #[cfg(any(debug_assertions, test))]
        {
            // k-mers repartition between dense and sparse index
            log::info!(
                "The index contains {:?} 'dense' k-mers and {:?} 'sparse' k-mers (total k-mers: {:?})",
                atomic_dense_kmers_count.get_mut().separate_with_commas(),
                atomic_sparse_kmers_count.get_mut().separate_with_commas(),
                total_kmers.get_mut().separate_with_commas()
            );

            let ones = atomic_sparse_one_seen.get_mut();
            let silent = *atomic_sparse_kmers_count.get_mut() - *ones;
            let fp = atomic_sparse_fp_seen.get_mut();
            // k-mers indexed in the sparse index, FP silent and FP seen
            log::info!(
                "Among the {:?} k-mers added in the 'sparse' index, {:?} encountered hash collisions ({:?} silent and {:?} misleading).",
                atomic_sparse_kmers_count.get_mut().separate_with_commas(),
                (silent + *fp).separate_with_commas(),
                silent.separate_with_commas(),
                fp.separate_with_commas()
            );
        }

        // Merge the chunk-based bloom filters, if needed
        if chunks.len() > 1 {
            merge::merge_all_partitions_of_chunks(
                &dir_path,
                &dir_path, // output
                partitioned_bf_size,
                parameters.abundance_number.get(),
                color_chunks,
                parameters.partition_number,
            )?;
        }

        // write metadata info to disk
        self.indexed_file_names = get_file_names(&file_paths);
        self.save_to_disk()?;

        Ok((file_paths, dir_path))
    }

    // TODO store bf_dir in Reindeer2 ?
    /// Performs a query on the index and writes the results on the disk.
    pub fn query(
        &self,
        fasta_file: &str,
        bf_dir: &str,
        output_file: &str,
        output_format: OutputFormat,
        coverage: f32,
    ) -> io::Result<()> {
        // let sampler = NoSampler::new();
        let sampler = &self.parameters.sampling_strategy;
        let reader = read_file(fasta_file)?;
        let mut writer = BufWriter::new(File::create(output_file)?);
        let output_format =
            query::EnrichedOutputFormat::from_pub_output_format(output_format, bf_dir);
        // write the header of the result file
        query::write_header(&self.indexed_file_names, &output_format, &mut writer)
            .expect("should have been able to write the header of the result file");
        let parameters = &self.parameters;

        let base = compute_base(parameters.abundance_number, parameters.abundance_max);

        // Process FASTA in chunks of 10 000 000 records
        process_fasta_in_batches(reader, 10_000_000, |batch| {
            let sequence_results = match sampler {
                None => self.query_single_fasta_batch(bf_dir, base, batch, &NoSampler::new()),
                Some(SamplingStrategy::KmerSampling { last_bits_to_zero }) => self
                    .query_single_fasta_batch(
                        bf_dir,
                        base,
                        batch,
                        &KmerSampler::new(*last_bits_to_zero),
                    ),
                Some(SamplingStrategy::MinimizerSampling { last_bits_to_zero }) => self
                    .query_single_fasta_batch(
                        bf_dir,
                        base,
                        batch,
                        &MinimizerSampler::new(*last_bits_to_zero),
                    ),
            };
            // Now `sequence_results` has the combined data for this batch.
            // Let's compute the output in the requested format.
            query::write_kmer_query(
                batch,
                &output_format,
                coverage,
                &sequence_results,
                &self.indexed_file_names,
                &mut writer,
            )
            .expect("should have been able to write the query result");
        })
        .expect("should have been able to process fasta files");

        Ok(())
    }

    // TODO use compile time to decide wether to sort the output or not
    fn query_single_fasta_batch<S>(
        &self,
        bf_dir: &str,
        base: f64,
        batch: &[fasta::Record],
        sampler: &S,
    ) -> Vec<Vec<Vec<ApproxAbundance>>>
    where
        S: Sampler,
    {
        let parameters = &self.parameters;
        let s = parameters.k - parameters.findere_z;
        let partition_smers = query::build_partitions_smers(
            batch,
            s,
            parameters.m,
            parameters.partition_number as u64,
            parameters.canonical,
            sampler,
        );

        // --- PARALLEL PHASE: process each partition's k-mers in parallel ---
        // This will store final results for *all* sequences in this batch.
        // Key: sequence header; Value: vector of vector of abundances
        let result_with_positions: HashMap<usize, Vec<Vec<(u32, ApproxAbundance)>>> =
            partition_smers
                .into_par_iter()
                // 1) Create a local HashMap in each thread
                .fold(
                    // OPTIMIZE use compile time monomorphization to get rid of this (usize, u16)
                    // (I expect to go from 128 per element in the vector to 16)
                    // (I needed this extra usize to get the correct order of kmer in the output)
                    HashMap::<usize, Vec<Vec<(u32, ApproxAbundance)>>>::new,
                    |local_results: HashMap<usize, Vec<Vec<(u32, ApproxAbundance)>>>,
                     (partition_index, smers)| {
                        query::fold_into_hashmap(
                            local_results,
                            partition_index,
                            smers,
                            base,
                            bf_dir,
                            parameters.bf_size,
                            parameters.partition_number,
                            parameters.nb_color,
                            parameters.abundance_number.get(),
                            parameters.dense_option,
                        )
                    },
                )
                // 2) Reduce all local HashMaps into a single HashMap
                .reduce(
                    HashMap::<usize, Vec<Vec<(u32, ApproxAbundance)>>>::new,
                    query::merge_results,
                );

        let mut smer_res: Vec<Vec<Vec<ApproxAbundance>>> = vec![vec![]; batch.len()];
        let empty = vec![vec![]; self.parameters.nb_color];
        (0..smer_res.len())
            .map(|header_id| {
                (
                    header_id,
                    result_with_positions.get(&header_id).unwrap_or(&empty),
                )
            })
            .for_each(|(header, color_vectors)| {
                smer_res[header] = color_vectors
                    .iter()
                    .map(|x| {
                        query::sort_abundance_vec(
                            x,
                            batch[header].seq().len() - s + 1,
                            #[cfg(any(debug_assertions, test))]
                            s,
                        )
                    })
                    .collect();
            });

        fimpera(smer_res, parameters.findere_z)
    }
}

fn get_file_names(file_paths: &[String]) -> Vec<String> {
    file_paths
        .iter()
        .map(|file_path| {
            // get the file name, excluding everything after the first "."
            let name = Path::new(file_path)
                .file_prefix()
                .unwrap_or_else(|| {
                    panic!(
                        "should have been able to extact the name of the file {}",
                        file_path
                    )
                })
                .to_str()
                .unwrap_or_else(|| {
                    panic!(
                        "should have been able to convert the filename {} to UTF-8",
                        file_path
                    )
                });
            String::from(name)
        })
        .collect()
}

/// Collects the iterator in batches and runs `process_batch` on each batch.
fn process_fasta_in_batches<R: io::BufRead>(
    reader: R,
    batch_size: usize,
    mut process_batch: impl FnMut(&[fasta::Record]),
) -> io::Result<()> {
    let fasta_reader = fasta::Reader::new(reader);
    let mut batch = Vec::with_capacity(batch_size);

    for result in fasta_reader.records() {
        let record = result?;
        batch.push(record);

        if batch.len() >= batch_size {
            process_batch(&batch);
            batch.clear();
        }
    }

    // final batch
    if !batch.is_empty() {
        process_batch(&batch);
    }

    Ok(())
}

// pub fn build_index_muset(
//     unitigs_file: String,
//     matrix_file: String,
//     kmer_size: usize,
//     minimizer_size: usize,
//     bf_size: u64,
//     partition_number: usize,
//     color_nb: usize,
//     abundance_number: usize,
//     abundance_max: u16,
//     output_dir: &str,
//     dense_option: bool,
//     threshold: usize,
//     canonical: bool,
//     debug: bool,
// ) -> io::Result<(Vec<String>, String)> {
//     let mut total_kmers = atomic::AtomicU64::new(0);
//     let mut atomic_dense_kmers_count = atomic::AtomicU64::new(0);
//     let mut atomic_sparse_kmers_count = atomic::AtomicU64::new(0);
//     let mut atomic_sparse_one_seen = atomic::AtomicU64::new(0);
//     let mut atomic_sparse_fp_seen = atomic::AtomicU64::new(0);
//     let kmer_counts_vector: Arc<Mutex<Vec<usize>>> = Arc::new(Mutex::new(vec![color_nb, 0]));
//     let base = compute_base(abundance_number, abundance_max);
//     let max_map_size = 1_000_000;
//     if debug {
//         println!("In debug mode... the tool may take (much) longer than usual.");
//         println!("Using log base {}", base);
//     }
//     println!("Initializing Bloom filter slices...");

//     let (_, dir_path) = create_dir_and_files(partition_number, output_dir)?;

//     // Shared data structures protected by Mutex for safe parallel access
//     let maybe_dense_indexes: Option<Arc<DenseIndex>> = if dense_option {
//         Some(Arc::new(DenseIndex::with_partition_number(
//             partition_number,
//         )))
//     } else {
//         None
//     };

//     let bloom_filters = Arc::new(Filters::with_number_partition(
//         partition_number,
//         color_nb,
//         bf_size as usize, // TODO check unit
//         abundance_number,
//     ));

//     let unitigs_reader = match read_file(&unitigs_file) {
//         Ok(r) => r,
//         Err(e) => {
//             panic!("Failed to open file {}: {}", unitigs_file, e);
//         }
//     };
//     let matrix_reader = match read_file(&matrix_file) {
//         Ok(r) => r,
//         Err(e) => {
//             panic!("Failed to open file {}: {}", matrix_file, e);
//         }
//     };
//     let len_matrix = matrix_reader.lines().count();
//     let len_unitigs = unitigs_reader.lines().count();
//     if len_matrix * 2 != len_unitigs {
//         panic!("The number of unitigs doesn't match the matrix size.");
//     }

//     let mut matrix_lines = read_file(&matrix_file)?.lines();
//     let mut unitigs_lines = read_file(&unitigs_file)?.lines();

//     let mut partition_kmers: HashMap<usize, Vec<(u64, u16, usize, usize)>> = HashMap::new(); // keep kmer info to fill BFs

//     let mut kmer_counts_vector_locked = kmer_counts_vector
//         .lock()
//         .expect("Failed to lock the kmer counter hashmap");
//     for _ in 0..len_matrix {
//         let abundance_line = matrix_lines.next().unwrap().unwrap();
//         let mut abundance_iter = abundance_line.trim().split(" ");
//         let unitig_id: &str = abundance_iter.next().unwrap();

//         let abundances = abundance_iter
//             .enumerate()
//             .map(|(i, ab_value_str)| {
//                 let ab_value = ab_value_str
//                     .parse::<f64>()
//                     .expect("Failed to parse the abundance values in the matrix")
//                     as u16;
//                 if ab_value == 0 {
//                     0
//                 } else {
//                     kmer_counts_vector_locked[i] += 1; // select the count corresponding to the right color
//                     (compute_log_abundance(ab_value, base, abundance_max) + 1) as u8
//                 }
//             })
//             .collect::<Vec<u8>>();

//         let tmp_abundance_vector = abundances.clone();
//         let number_of_zeros: usize = tmp_abundance_vector.into_iter().filter(|x| *x == 0).count();

//         let unitig_line = unitigs_lines.next().unwrap().unwrap();
//         if &unitig_line.split(" ").next().unwrap()[1..] != unitig_id {
//             println!(
//                 "{} vs {}",
//                 &unitig_line.split(" ").next().unwrap()[1..],
//                 unitig_id
//             );
//             panic!("The unitig id dooesn't match between the unitigs file and the matrix file.")
//         }

//         let unitig_line = unitigs_lines.next().unwrap().unwrap();
//         let seq_str = unitig_line.trim();
//         for (kmer_hash, minimizer) in
//             kmer_minimizers_seq_level(seq_str.as_bytes(), kmer_size, minimizer_size, canonical)
//         {
//             let partition = (minimizer % (partition_number as u64)) as usize;
//             let new_kmers_added = (abundances.len() - number_of_zeros) as u64;
//             total_kmers.fetch_add(new_kmers_added, atomic::Ordering::Relaxed);

//             if dense_option && number_of_zeros <= threshold {
//                 match &maybe_dense_indexes {
//                     Some(dense_indexes) => {
//                         dense_indexes.insert_abundance_to_partition(
//                             partition,
//                             kmer_hash,
//                             abundances.clone(),
//                         );
//                         atomic_dense_kmers_count
//                             .fetch_add(new_kmers_added, atomic::Ordering::Relaxed);
//                     }
//                     None => panic!("Failed to build the dense index."),
//                 }
//             } else {
//                 atomic_sparse_kmers_count.fetch_add(new_kmers_added, atomic::Ordering::Relaxed);
//                 let kmers_entry = partition_kmers // separate the kmers per partition
//                     .entry(partition)
//                     .or_default();

//                 for (path_num, log_plusone) in abundances.iter().enumerate() {
//                     let real_log_abundance = (log_plusone - 1) as u16;
//                     if real_log_abundance > 0 {
//                         kmers_entry.push((
//                             kmer_hash,
//                             real_log_abundance,
//                             path_num,
//                             0, // chunk index
//                         ))
//                     }
//                 }
//             }
//             if partition_kmers.len() >= max_map_size {
//                 bloom_filters.extend_by_draining_partitions_map(&mut partition_kmers);
//             }
//         }
//     }
//     // Index remaining kmers
//     bloom_filters.extend_by_draining_partitions_map(&mut partition_kmers);

//     #[cfg(any(debug_assertions, test))]
//     bloom_filters.update_sparse_counts(
//         &atomic_sparse_one_seen,
//         &atomic_sparse_fp_seen,
//         abundance_number,
//     )?;

//     if let Err(e) = bloom_filters.write_to_disk(&dir_path, &vec![color_nb], partition_number, 0) {
//         eprintln!("Error writing Bloom filters on disk: {}", e);
//     }

//     // After processing all chunks, write the dense indexes to disk
//     if let Some(dense_indexes) = maybe_dense_indexes {
//         dense_indexes.write_to_disk(&dir_path)?;
//     }

//     write_kmer_counts_to_disk(&dir_path, &kmer_counts_vector)?;

//     if debug {
//         // k-mers repartition between dense and sparse index
//         println!(
//             "The index contains {:?} 'dense' k-mers and {:?} 'sparse' k-mers (total k-mers: {:?})",
//             atomic_dense_kmers_count.get_mut().separate_with_commas(),
//             atomic_sparse_kmers_count.get_mut().separate_with_commas(),
//             total_kmers.get_mut().separate_with_commas()
//         );

//         let ones = atomic_sparse_one_seen.get_mut();
//         let silent = *atomic_sparse_kmers_count.get_mut() - *ones;
//         let fp = atomic_sparse_fp_seen.get_mut();
//         // k-mers indexed in the sparse index, FP silent and FP seen
//         println!("Among the {:?} k-mers added in the 'sparse' index, {:?} encountered hash collisions ({:?} silent and {:?} misleading).",
//             atomic_sparse_kmers_count.get_mut().separate_with_commas(),
//             (silent+*fp).separate_with_commas(),
//             silent.separate_with_commas(),
//             fp.separate_with_commas());
//     }

//     // Merge the chunk-based bloom filters, if needed
//     // If there was only one chunk, rename files directly
//     for partition_idx in 0..partition_number {
//         let input_path = format!(
//             "{}/partition_bloom_filters_c0_p{}.bin",
//             dir_path, partition_idx
//         );
//         let output_path = format!(
//             "{}/partition_bloom_filters_p{}.bin",
//             dir_path, partition_idx
//         );
//         std::fs::rename(&input_path, &output_path)?;
//     }

//     // write partition info to a CSV or your desired format
//     let _ = write_partition_to_csv(
//         &dir_path,
//         kmer_size,
//         minimizer_size,
//         bf_size,
//         partition_number,
//         color_nb,
//         abundance_number,
//         abundance_max,
//         dense_option,
//         canonical,
//     );

//     Ok((vec!["".to_string()], dir_path))
// }

// === GENERAL ===

// --- BF MANAGEMENT ---

const fn compute_base_position(
    smer_hash: u64,
    partitioned_bf_size: usize,
    color_number: usize,
    abundance_number: usize,
) -> u64 {
    // return the first position of the concerned column in the partition
    let position = smer_hash % (partitioned_bf_size as u64);
    position * (color_number as u64) * (abundance_number as u64)
}

// fn _get_current_chunk_index(i: usize, chunk_sizes: &Vec<usize>, partition_nb: usize) -> usize {
//     let mut cumulative_size = 0;

//     for (chunk_idx, &_chunk_size) in chunk_sizes.iter().enumerate() {
//         cumulative_size += partition_nb;
//         if i < cumulative_size {
//             return chunk_idx;
//         }
//     }
//     panic!(
//         "Index {} out of bounds for chunk sizes {:?}",
//         i, chunk_sizes
//     );
// }

// // this part fills the BFs per partition
// fn flush_map_into_bfs(
//     partition_kmers: &mut HashMap<usize, Vec<(u64, u16, usize, usize)>>,
//     bloom_filters: &Arc<Filters>,
// ) {
//     // iterates and empties the hash map when needed
//     for (partition_index, kmers) in partition_kmers.drain() {
//         bloom_filters.extend_partition(partition_index, &kmers);
//     }
// }

// fn flush_map_into_bfs(
//     partition_kmers: &mut HashMap<usize, Vec<(u64, u16, usize, usize)>>,
//     bloom_filters: &Arc<Vec<Mutex<RoaringBitmap>>>,
//     partitioned_bf_size: usize,
//     color_number: usize,
//     abundance_number: usize,
// ) -> io::Result<()> {
//     // this part fills the BFs per paritition
//     for (partition_index, kmers) in partition_kmers.drain() {
//         // iterates and empties the hash map when needed
//         let mut kmer_hashes_to_update = Vec::new();
//         for (kmer_hash, log_abundance, path_num, _chunk_index) in kmers {
//             // select the bit in the BF
//             let position = compute_location_filter(
//                 kmer_hash,
//                 partitioned_bf_size,
//                 color_number,
//                 path_num,
//                 abundance_number,
//                 log_abundance,
//             );
//             kmer_hashes_to_update.push(position as u32); // accumulate bits to be modified for this bf
//         }
//         let mut bloom_filter = bloom_filters[partition_index] // select the correct BF for the given partition
//             .lock()
//             .expect("Failed to lock bloom filter");
//         bloom_filter.extend(kmer_hashes_to_update);
//     }
//     Ok(())
// }

// --- WRITE INDEX ON DISK ---

fn create_dir_and_files(
    num_partition: usize,
    output_dir: &str,
) -> io::Result<(Vec<String>, String)> {
    log::info!("Writing partitioned files in directory: {}", output_dir);
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

fn write_kmer_counts_to_disk(
    dir_path: &str,
    kmer_counts_vector: &Arc<Mutex<Vec<usize>>>,
) -> io::Result<()> {
    // OPTIMIZE we may be able to drop the lock before writing to disk
    let file_path = Path::new(dir_path).join("kmer_counts_per_color.bin");
    let file = File::create(&file_path)?;
    let mut writer = BufWriter::new(file);
    let mut locked_vector = kmer_counts_vector.lock().expect(
        "fatal error: a thread holding the mutex panicked, so this thread will panic as well",
    );
    let locked_vector_ref: &Vec<usize> = &locked_vector;
    let binary_encoded = bincode::serialize(locked_vector_ref)
        .expect("should have been able to serialize the count of k-mers");
    locked_vector.clear();
    drop(locked_vector);
    writer.write_all(&binary_encoded)?;
    Ok(())
}

// --- LOAD FROM DISK ---

// /// Rload index metadata from the CSV.
// // TODO why output_csv ?
// // TODO we could return a Reindeer2 directly
// // TODO why is 0 the default ?
// fn read_partition_from_csv(
//     bf_dir: &str,
//     output_csv: &str,
// ) -> io::Result<(usize, usize, u64, usize, usize, usize, u16, bool, bool)> {
//     let csv_path = format!("{}/{}", bf_dir, output_csv);
//     let mut reader = csv::Reader::from_reader(BufReader::new(File::open(csv_path)?));
//     let values = reader.records().next().expect("Index CSV is empty")?;
//     let k = values[0].parse().unwrap_or(0);
//     let m = values[1].parse().unwrap_or(0);
//     let bf_size = values[2].parse().unwrap_or(0);
//     let partition_number = values[3].parse().unwrap_or(0);
//     let color_number = values[4].parse().unwrap_or(0);
//     let abundance_number = values[5].parse().unwrap_or(0);
//     let abundance_max = values[6].parse().unwrap_or(0);
//     let dense_option = values[7].parse().unwrap_or(false);
//     let canonical = values[8].parse().unwrap_or(true);

//     Ok((
//         k,
//         m,
//         bf_size,
//         partition_number,
//         color_number,
//         abundance_number,
//         abundance_max,
//         dense_option,
//         canonical,
//     ))
// }

// --- FOF MANAGEMENT ---

/// Reads a file of file and return each lines, along with the number of lines in it.
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

fn split_fof(lines: &[String], number_of_chunks: usize) -> io::Result<Vec<Vec<String>>> {
    let number_colors = lines.len();

    // Determine split_factor
    let split_factor = if number_colors < number_of_chunks {
        1
    } else {
        number_colors.div_ceil(number_of_chunks)
    };

    let mut fof_chunks = vec![vec![]; split_factor];

    for (i, line) in lines.iter().enumerate() {
        let chunk_index = i / number_of_chunks;
        fof_chunks[chunk_index].push(line.clone());
    }

    Ok(fof_chunks)
}

// // TODO discuss: used to return a Result, but only the OK variant was returnd
// pub fn explore_muset_dir(dir_str: &str) -> (String, String, usize) {
//     let dir_path = Path::new(dir_str);
//     if !dir_path.is_dir() {
//         panic!("{} is not a directory", dir_str);
//     }
//     let unitigs_path = dir_path.join("unitigs.fa");
//     if !unitigs_path.exists() {
//         panic!("File not found : {:#?}", unitigs_path);
//     }
//     let abundance_path = dir_path.join("unitigs.abundance.mat");
//     if !abundance_path.exists() {
//         panic!("File not found : {:#?}", abundance_path);
//     }

//     let reader = BufReader::new(
//         File::open(&abundance_path).expect(&format!("Failed to open {:#?}", abundance_path)),
//     );
//     let err = &format!("Failed to read {:#?}", abundance_path);
//     let firstline = reader.lines().next().expect(err).expect(err);

//     // number of color == number of \t in the first line ?
//     // TODO check it's the same
//     let color_nb = firstline.chars().filter(|c| *c == '\t').count();
//     (
//         unitigs_path.to_string_lossy().to_string(),
//         abundance_path.to_string_lossy().to_string(),
//         color_nb,
//     )
// }

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
        Err(io::Error::other(format!(
            "Header {} does not contain a recognized count field (km:f: or ka:f:)",
            header,
        )))
    }
}

fn extract_count(header: &str, header_type: &HeaderType) -> Result<u16, io::Error> {
    match header_type {
        HeaderType::BCalm => extract_count_from_bcalm_header(header),
        HeaderType::Logan => extract_count_from_logan_header(header),
    }
}

fn extract_count_from_bcalm_header(header: &str) -> Result<u16, io::Error> {
    header
        .split_whitespace()
        .find(|&part| part.starts_with("km:f:"))
        .and_then(|km_part| km_part.trim_start_matches("km:f:").parse::<f32>().ok())
        .map(|float_val| float_val as u16) //any value over max u16 will be clamped to the maximum possible value
        .ok_or_else(|| io::Error::other("km:f: value not found in header or could not be parsed"))
}

fn extract_count_from_logan_header(header: &str) -> Result<u16, io::Error> {
    header
        .split_whitespace()
        .find(|&part| part.starts_with("ka:f:"))
        .and_then(|ka_part| ka_part.trim_start_matches("ka:f:").parse::<f32>().ok())
        .map(|float_val| float_val as u16)
        .ok_or_else(|| io::Error::other("ka:f: value not found in header or could not be parsed"))
}

// --- ABUNDANCE ENCODING ---

fn compute_log_abundance(value: NonZero<u16>, base: f64, max: NonZero<u16>) -> u16 {
    assert!(base > 0.0, "base must be greater than 0");
    let value = value.get();
    let max = max.get();

    let value_f = value.min(max) as f64;
    let threshold = 1.0 / (base - 1.0);

    if value_f < threshold {
        value - 1
    } else {
        (value_f.ln() / base.ln() + (base - 1.0).ln() / base.ln() + threshold - 1.0) as u16
    }
}

fn approximate_value(log_value: u16, base: f64) -> u16 {
    if base <= 0.0 {
        panic!("base must be greater than 0");
    }
    let threshold = 1.0 / (base - 1.0);
    let logf = log_value as f64;
    if logf < threshold {
        log_value + 1
    } else {
        (base.powf((logf + 1.0) - threshold) * threshold) as u16
    }
}

// TOUN
fn compute_base(abundance_number: NonZero<usize>, abundance_max: NonZero<u16>) -> f64 {
    let abundance_numberf = abundance_number.get() as f64;
    const TOL: f64 = 1e-9;
    let abundance_maxf = abundance_max.get() as f64;

    let equation = |b: f64| -> f64 {
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

#[cfg(test)]
mod tests {

    use crate::reindeer2::{merge::merge_multiple_indexes, test_utils::AutoRemoveDirectory};

    use super::*;

    use itertools::Itertools;
    use rstest::{fixture, rstest};
    use rstest_reuse::{self, *};

    #[test]
    fn test_read_fof_file() {
        let test_file_path = "test_fof.txt";

        // Create the test file
        {
            let mut test_file =
                File::create(test_file_path).expect("Failed to create test FOF file");
            writeln!(test_file, "/home/test/path/file1.fasta")
                .expect("Failed to write to test FOF file");
            writeln!(test_file, "/home/test/path/file2.fasta")
                .expect("Failed to write to test FOF file");
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
        assert_eq!(result, expected);
    }
    #[test]
    fn test_extract_count_from_logan_header4() {
        let header = ">0 ka:f:65537   L:+:5:+ L:+:5392806:+  L:-:1:-";
        let expected = 65535;
        let result = extract_count_from_logan_header(header).unwrap();
        assert_eq!(result, expected);
    }

    #[test]
    fn test_compute_log_abundance_gives_zero() {
        let value = NonZero::new(1).unwrap();
        let max = NonZero::new(65535).unwrap();

        let result = compute_log_abundance(value, 2.0, max);
        assert_eq!(result, 0, "expected log for 1 to be 0");
        let result2 = compute_log_abundance(value, 1.5, max);
        assert_eq!(result2, 0, "expected log for 1 to be 0");
    }

    #[test]
    fn test_compute_log_abundance_is_linear() {
        let max = NonZero::new(65535).unwrap();
        let b = 1.1;
        for n in 1..9 {
            let n = NonZero::new(n).unwrap();
            let result = compute_log_abundance(n, b, max);
            assert_eq!(
                result,
                n.get() - 1,
                "expected abundance function to be linear for the first values"
            );
        }
        assert_eq!(
            compute_log_abundance(NonZero::new(10).unwrap(), b, max),
            compute_log_abundance(NonZero::new(11).unwrap(), b, max),
            "expected abundance function not to be linear at the threshold"
        );
    }

    #[test]
    fn test_compute_log_abundance_positive_values() {
        let value = NonZero::new(8).unwrap();
        let max = NonZero::new(65535).unwrap();
        let result = compute_log_abundance(value, 2.0, max);
        assert!(
            ((result as f64) - 3.0).abs() < 1e-6,
            "expected log base 2 of 8 to be 3"
        );
    }

    #[test]
    fn test_compute_log_abundance_with_large_values() {
        let value = NonZero::new(65535).unwrap();
        let max = value;
        let result = compute_log_abundance(value, 2.0, max);
        let expected = (value.get() as f64).log2().floor() as u16;
        assert_eq!(
            result, expected,
            "expected log base 2 of 65535 to be {}, but got {}",
            expected, result
        );
    }

    #[test]
    fn test_compute_log_abundance_above_max() {
        let value = NonZero::new(65535).unwrap();
        let max = NonZero::new(256).unwrap();
        let result = compute_log_abundance(value, 2.0, max);
        let expected = 256_f64.log2().floor() as u16;
        assert_eq!(
            result, expected,
            "expected abundance to be scaled at log base 2 of 256 should have been {}, but got {}",
            expected, result
        );
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
    #[should_panic(expected = "base must be greater than 0")]
    fn test_compute_log_abundance_non_positive_base() {
        let value = NonZero::new(8).unwrap();
        let max = NonZero::new(65535).unwrap();
        compute_log_abundance(value, 0.0, max);
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
        assert_eq!(
            result2 / 2,
            result1,
            "check approximation consistency with base sqrt(2)"
        );
        assert_eq!(
            result3 / 2,
            result2,
            "check approximation consistency with base sqrt(2)"
        );
    }

    #[test]
    fn test_compute_base_with_positive_values() {
        let abundance_number = NonZero::new(16).unwrap();
        let abundance_max = NonZero::new(1024).unwrap();
        let result = compute_base(abundance_number, abundance_max);
        assert!(
            (result - 1.5635206).abs() < 1e-6,
            "expected base for 16 with max 1024 to be ~1.56"
        );
    }

    #[test]
    fn test_compute_base_with_large_abundance_number() {
        let abundance_number = NonZero::new(32).unwrap();
        let abundance_max = NonZero::new(1024).unwrap();
        let result = compute_base(abundance_number, abundance_max);
        assert!(
            (result - 1.218096).abs() < 1e-6,
            "expected base for 32 with max 1024 to be approximately ~1.22"
        );
    }

    #[test]
    fn test_compute_base_with_small_abundance_number() {
        let abundance_number = NonZero::new(8).unwrap();
        let abundance_max = NonZero::new(1024).unwrap();
        let result = compute_base(abundance_number, abundance_max);
        assert!(
            (result - 2.740397).abs() < 1e-6,
            "expected base for 8 with max 1024 to be ~2.74"
        );
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

    #[template]
    #[rstest]
    #[case(0)]
    #[case(2)]
    #[case(4)]
    #[case(6)]
    #[case(8)]
    #[case(10)]
    fn findere_fixture(#[case] z: usize) {}

    #[fixture]
    pub fn random_directory() -> AutoRemoveDirectory {
        AutoRemoveDirectory::create_random()
    }

    pub fn wrong_capacity() -> usize {
        // this is used to give a wrong capacity, but it should no matter for tests
        // TODO change it later
        0
    }

    #[apply(findere_fixture)]
    fn test_write_and_read_metadata(#[case] z: usize, random_directory: AutoRemoveDirectory) {
        let bf_dir = random_directory.filename().to_str().unwrap();
        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024,
            partition_number: 4,
            nb_color: 3,
            abundance_number: NonZero::new(2).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(512).unwrap(),
            dense_option: false,
            canonical: false,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let indexed_file_names = vec![String::from("a"), String::from("b")];

        let mut expected = Reindeer2::new(parameters, String::from(bf_dir));
        expected.set_indexed_file_names(indexed_file_names);
        fs::create_dir_all(bf_dir).expect("Failed to create test directory");
        expected.save_to_disk().unwrap();

        let actual = Reindeer2::load_from_disk(bf_dir).unwrap();
        assert_eq!(actual, expected);
    }

    /// Helper function for tests
    /// Creates an index indexing `file_paths`, load it from disk, and query file `file_paths[query_file_id]`.
    fn create_build_query(
        parameters: Parameters,
        chunks_size: usize,
        threshold: usize,
        file_paths: Vec<String>,
        query_file_id: usize,
        test_dir: impl Into<String>,
    ) -> String {
        let mut index = Reindeer2::new(parameters, test_dir.into());
        let (_file_paths, index_dir) = index
            .build(file_paths.clone(), chunks_size, threshold, false)
            .expect("Failed to build index");

        let query_results_path = format!("{}/query_results.csv", index_dir);

        let index_from_disk =
            Reindeer2::load_from_disk(&index_dir).expect("Failed to load index infos from disk");
        index_from_disk
            .query(
                &file_paths[query_file_id],
                &index_dir,
                &query_results_path,
                OutputFormat::Median { normalized: None },
                0.5,
            )
            .expect("Failed to query sequences");
        query_results_path
    }

    fn load_query_result_csv<P: AsRef<Path>>(path: P) -> Vec<(String, String, usize)> {
        let file = File::open(&path).expect("Failed to open query results");

        csv::Reader::from_reader(file)
            .records()
            .map(|record| record.expect("Failed to read record"))
            .map(|record| {
                let header = record[0].to_string(); // Now reading the header instead of the sequence ID
                let color = record[1].parse().expect("Failed to parse color's name");
                let abundance = record[2].parse().expect("Failed to parse abundance");
                (header, color, abundance)
            })
            .collect()
    }

    #[apply(findere_fixture)]
    fn test_build_and_query_index_single_sparse(
        #[case] z: usize,
        random_directory: AutoRemoveDirectory,
    ) {
        let test_dir = random_directory.filename().to_str().unwrap();
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

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 1,
            abundance_number: NonZero::new(256).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: false,
            canonical: false,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let mut index = Reindeer2::new(parameters, String::from(test_dir));
        let (_file_paths, index_dir) = index
            .build(vec![file1_path.clone()], chunks_size, threshold, false)
            .expect("Failed to build index");

        let query_results_path = format!("{}/query_results.csv", index_dir);

        let index_from_disk =
            Reindeer2::load_from_disk(&index_dir).expect("Failed to load index infos from disk");
        index_from_disk
            .query(
                &file1_path,
                &index_dir,
                &query_results_path,
                OutputFormat::Median { normalized: None },
                0.5,
            )
            .expect("Failed to query sequences");

        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [
            (">seq1 ka:f:30".to_string(), String::from("file1Q"), 29),
            (">seq2 ka:f:30".to_string(), String::from("file1Q"), 29),
            (">seq3 ka:f:2".to_string(), String::from("file1Q"), 2), // Values with errors due to log conversion
        ];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
    }

    #[apply(findere_fixture)]
    fn test_build_and_query_index_single_dense(
        #[case] z: usize,
        random_directory: AutoRemoveDirectory,
    ) {
        let test_dir = random_directory.filename().to_str().unwrap();
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

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 1,
            abundance_number: NonZero::new(255).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: true,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let query_results_path = create_build_query(
            parameters,
            chunks_size,
            threshold,
            vec![file1_path],
            0,
            test_dir,
        );

        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [
            (">seq1 ka:f:30".to_string(), String::from("file1Q"), 29),
            (">seq2 ka:f:30".to_string(), String::from("file1Q"), 29),
            (">seq3 ka:f:2".to_string(), String::from("file1Q"), 2), // Values with errors due to log conversion
        ];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
    }

    #[apply(findere_fixture)]
    fn test_build_and_query_index_sparse(#[case] z: usize, random_directory: AutoRemoveDirectory) {
        let test_dir = random_directory.filename().to_str().unwrap();
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

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 2,
            abundance_number: NonZero::new(256).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: false,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let query_results_path = create_build_query(
            parameters,
            chunks_size,
            threshold,
            vec![file1_path, file2_path],
            0,
            test_dir,
        );
        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [
            (">seq1 ka:f:30".to_string(), String::from("file1Q"), 29),
            (">seq2 ka:f:30".to_string(), String::from("file1Q"), 29),
            (">seq3 ka:f:2".to_string(), String::from("file1Q"), 2),
            (">seq3 ka:f:2".to_string(), String::from("file2Q"), 997),
        ];
        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
    }

    #[apply(findere_fixture)]
    fn test_build_and_query_index_dense(#[case] z: usize, random_directory: AutoRemoveDirectory) {
        let test_dir = random_directory.filename().to_str().unwrap();
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

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 2,
            abundance_number: NonZero::new(256).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: true,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let query_results_path = create_build_query(
            parameters,
            chunks_size,
            threshold,
            vec![file1_path, file2_path],
            0,
            test_dir,
        );

        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [
            (">seq1 ka:f:30".to_string(), String::from("file1Q"), 29),
            (">seq2 ka:f:30".to_string(), String::from("file1Q"), 29),
            (">seq3 ka:f:2".to_string(), String::from("file1Q"), 2),
            (">seq3 ka:f:2".to_string(), String::from("file2Q"), 997),
        ];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
    }

    #[apply(findere_fixture)]
    fn test_build_and_query_index_longseq_simple_sparse(
        #[case] z: usize,
        random_directory: AutoRemoveDirectory,
    ) {
        let test_dir = random_directory.filename().to_str().unwrap();
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
            writeln!(file2, "AAAAATGATAGTAGAAAAAAATTTTAAAAAAACACCCCTGG")
                .expect("Failed to write sequence");
        }

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 2,
            abundance_number: NonZero::new(256).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: false,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let query_results_path = create_build_query(
            parameters,
            chunks_size,
            threshold,
            vec![file1_path, file2_path],
            1,
            test_dir,
        );

        // Validate the results written to the query results CSV file
        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [(">seq3 ka:f:1000".to_string(), String::from("file2Q"), 997)];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
    }

    #[apply(findere_fixture)]
    fn test_build_and_query_index_longseq_simple_dense(
        #[case] z: usize,
        random_directory: AutoRemoveDirectory,
    ) {
        let test_dir = random_directory.filename().to_str().unwrap();
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
            writeln!(file2, "AAAAATGATAGTAGAAAAAAATTTTAAAAAAACACCCCTGG")
                .expect("Failed to write sequence");
        }

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 2,
            abundance_number: NonZero::new(256).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: true,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let query_results_path = create_build_query(
            parameters,
            chunks_size,
            threshold,
            vec![file1_path, file2_path],
            1,
            test_dir,
        );

        // Validate the results written to the query results CSV file
        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [(">seq3 ka:f:1000".to_string(), String::from("file2Q"), 997)];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
    }

    #[apply(findere_fixture)]
    fn test_build_and_query_index2_sparse(#[case] z: usize, random_directory: AutoRemoveDirectory) {
        use std::io::Write;

        let test_dir = random_directory.filename().to_str().unwrap();
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

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 2,
            abundance_number: NonZero::new(255).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: false,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = 1;

        let query_results_path = create_build_query(
            parameters,
            chunks_size,
            threshold,
            vec![file1_path, file2_path],
            1,
            test_dir,
        );

        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [
            (">seq4 ka:f:1000".to_string(), String::from("file2Q"), 979),
            (">seq5 ka:f:1000".to_string(), String::from("file2Q"), 979),
        ];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
    }

    #[apply(findere_fixture)]
    fn test_build_and_query_index2_dense(#[case] z: usize, random_directory: AutoRemoveDirectory) {
        use std::io::Write;

        let test_dir = random_directory.filename().to_str().unwrap();
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

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 2,
            abundance_number: NonZero::new(255).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: true,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = 1;

        let query_results_path = create_build_query(
            parameters,
            chunks_size,
            threshold,
            vec![file1_path, file2_path],
            1,
            test_dir,
        );

        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [
            (">seq4 ka:f:1000".to_string(), String::from("file2Q"), 979),
            (">seq5 ka:f:1000".to_string(), String::from("file2Q"), 979),
        ];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
    }

    #[apply(findere_fixture)]
    fn test_build_and_query_index_sharedk_sparse(
        #[case] z: usize,
        random_directory: AutoRemoveDirectory,
    ) {
        let test_dir = random_directory.filename().to_str().unwrap();
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

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 2,
            abundance_number: NonZero::new(255).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(255).unwrap(),
            dense_option: false,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let query_results_path = create_build_query(
            parameters,
            chunks_size,
            threshold,
            vec![file1_path, file2_path],
            0,
            test_dir,
        );

        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [
            (">seq1 ka:f:30".to_string(), String::from("file1Q"), 30),
            (">seq2 ka:f:30".to_string(), String::from("file1Q"), 30),
            (">seq3 ka:f:1500".to_string(), String::from("file1Q"), 255), //because of abundance_max
            (">seq3 ka:f:1500".to_string(), String::from("file2Q"), 4),
        ];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
    }

    #[apply(findere_fixture)]
    fn test_build_and_query_index_sharedk_dense(
        #[case] z: usize,
        random_directory: AutoRemoveDirectory,
    ) {
        let test_dir = random_directory.filename().to_str().unwrap();
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

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 2,
            abundance_number: NonZero::new(255).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(255).unwrap(),
            dense_option: true,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let query_results_path = create_build_query(
            parameters,
            chunks_size,
            threshold,
            vec![file1_path, file2_path],
            0,
            test_dir,
        );

        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [
            (">seq1 ka:f:30".to_string(), String::from("file1Q"), 30),
            (">seq2 ka:f:30".to_string(), String::from("file1Q"), 30),
            (">seq3 ka:f:1500".to_string(), String::from("file1Q"), 255), //because of abundance_max
            (">seq3 ka:f:1500".to_string(), String::from("file2Q"), 4),
        ];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
    }

    #[apply(findere_fixture)]
    fn test_build_and_query_index_sharedk_sparse_smallchunks(
        #[case] z: usize,
        random_directory: AutoRemoveDirectory,
    ) {
        let test_dir = random_directory.filename().to_str().unwrap();
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

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 2,
            abundance_number: NonZero::new(255).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(255).unwrap(),
            dense_option: false,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 1;
        let threshold = parameters.nb_color;

        let query_results_path = create_build_query(
            parameters,
            chunks_size,
            threshold,
            vec![file1_path, file2_path],
            0,
            test_dir,
        );

        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [
            (">seq1 ka:f:30".to_string(), String::from("file1Q"), 30),
            (">seq2 ka:f:30".to_string(), String::from("file1Q"), 30),
            (">seq3 ka:f:1500".to_string(), String::from("file1Q"), 255), //because of abundance_max
            (">seq3 ka:f:1500".to_string(), String::from("file2Q"), 4),
        ];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
    }

    #[apply(findere_fixture)]
    fn test_build_and_query_index_sharedk_dense_smallchunks(
        #[case] z: usize,
        random_directory: AutoRemoveDirectory,
    ) {
        let test_dir = random_directory.filename().to_str().unwrap();
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

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 2,
            abundance_number: NonZero::new(255).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(255).unwrap(),
            dense_option: true,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 1;
        let threshold = parameters.nb_color;

        let query_results_path = create_build_query(
            parameters,
            chunks_size,
            threshold,
            vec![file1_path, file2_path],
            0,
            test_dir,
        );

        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [
            (">seq1 ka:f:30".to_string(), String::from("file1Q"), 30),
            (">seq2 ka:f:30".to_string(), String::from("file1Q"), 30),
            (">seq3 ka:f:1500".to_string(), String::from("file1Q"), 255), //because of abundance_max
            (">seq3 ka:f:1500".to_string(), String::from("file2Q"), 4),
        ];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
    }

    #[apply(findere_fixture)]
    fn test_build_and_query_index_long_sparse(
        #[case] z: usize,
        random_directory: AutoRemoveDirectory,
    ) {
        let test_dir = random_directory.filename().to_str().unwrap();
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
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACACAGATCAGAG")
                .expect("Failed to write sequence");
            writeln!(file1, ">seq2 ka:f:8").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT").expect("Failed to write sequence");
        }

        {
            let mut file2 = File::create(&file2_path).expect("Failed to create file2.fasta");
            writeln!(file2, ">seq4 ka:f:1000").expect("Failed to write header");
            writeln!(file2, "CAAAAAAAAAAAAAAAAAAAAACACCCCTGGAC").expect("Failed to write sequence");
        }

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 2,
            abundance_number: NonZero::new(255).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: false,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let query_results_path = create_build_query(
            parameters,
            chunks_size,
            threshold,
            vec![file1_path, file2_path],
            0,
            test_dir,
        );

        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [
            (">seq1 ka:f:30".to_string(), String::from("file1Q"), 29),
            (">seq2 ka:f:8".to_string(), String::from("file1Q"), 8),
        ];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
    }

    #[apply(findere_fixture)]
    fn test_build_and_query_index_long_dense(
        #[case] z: usize,
        random_directory: AutoRemoveDirectory,
    ) {
        let test_dir = random_directory.filename().to_str().unwrap();
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
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAAACACAGATCAGAG")
                .expect("Failed to write sequence");
            writeln!(file1, ">seq2 ka:f:8").expect("Failed to write header");
            writeln!(file1, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT").expect("Failed to write sequence");
        }

        {
            let mut file2 = File::create(&file2_path).expect("Failed to create file2.fasta");
            writeln!(file2, ">seq4 ka:f:1000").expect("Failed to write header");
            writeln!(file2, "CAAAAAAAAAAAAAAAAAAAAACACCCCTGGAC").expect("Failed to write sequence");
        }

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 2,
            abundance_number: NonZero::new(255).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: true,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let query_results_path = create_build_query(
            parameters,
            chunks_size,
            threshold,
            vec![file1_path, file2_path],
            0,
            test_dir,
        );

        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [
            (">seq1 ka:f:30".to_string(), String::from("file1Q"), 29),
            (">seq2 ka:f:8".to_string(), String::from("file1Q"), 8),
        ];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
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
        let base = 2.0;

        let mut bloom_filters: Vec<RoaringBitmap> = vec![RoaringBitmap::new(); partition_number];

        let partition_index_insert = (minimizer % (partition_number as u64)) as usize;
        let position_to_write = Filters::compute_location(
            kmer_hash,
            partitioned_bf_size,
            color_number,
            path_color_number,
            abundance_number,
            log_abundance,
        );
        bloom_filters[partition_index_insert].insert(position_to_write);

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
                position_to_write as u64,
                base_position + (log_abundance as u64),
                "Position mismatch: insert = {}, query = {}",
                position_to_write,
                base_position + (log_abundance as u64)
            );

            query::update_color_abundances(
                bitmap,
                base_position,
                color_number,
                abundance_number,
                0,
                base,
                &mut color_abundances,
            );
        }

        let results = color_abundances
            .into_iter()
            .enumerate()
            .filter_map(|(color, abundances)| {
                let min_abundance = abundances
                    .into_iter()
                    .map(|(_kmer_pos, abundance)| abundance)
                    .min();
                min_abundance.map(|abundance| (color, abundance))
            })
            .collect_vec();

        let expected_results = vec![
            (0, ApproxAbundance::from_dense(log_abundance as u8, base)),
            (1, ApproxAbundance::from_dense(0, base)),
            (2, ApproxAbundance::from_dense(0, base)),
        ]; // expect color 0 with abundance level 1

        assert_eq!(results, expected_results);
    }

    #[rstest]
    fn test_split_fof(random_directory: AutoRemoveDirectory) -> std::io::Result<()> {
        use std::fs::{self, File};
        use std::io::{BufRead, BufReader, Write};

        let test_dir = random_directory.filename().to_str().unwrap();
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

        // Test with chunks_size = 2
        let chunks_size = 2;
        let result = split_fof(&lines, chunks_size);

        assert!(result.is_ok(), "split_fof returned an error");

        let fof_chunks = result.unwrap();

        let expected_chunks = vec![
            vec![
                "/path/to/file1.fasta".to_string(),
                "/path/to/file2.fasta".to_string(),
            ],
            vec![
                "/path/to/file3.fasta".to_string(),
                "/path/to/file4.fasta".to_string(),
            ],
            vec!["/path/to/file5.fasta".to_string()],
        ];

        assert_eq!(
            fof_chunks, expected_chunks,
            "Mismatch in chunk distribution: expected {:?}, got {:?}",
            expected_chunks, fof_chunks
        );

        // Test with chunks_size = 5
        let chunks_size = 5;
        let result = split_fof(&lines, chunks_size);

        assert!(result.is_ok(), "split_fof returned an error");

        let fof_chunks = result.unwrap();

        let expected_chunks = vec![vec![
            "/path/to/file1.fasta".to_string(),
            "/path/to/file2.fasta".to_string(),
            "/path/to/file3.fasta".to_string(),
            "/path/to/file4.fasta".to_string(),
            "/path/to/file5.fasta".to_string(),
        ]];

        assert_eq!(
            fof_chunks, expected_chunks,
            "Mismatch in chunk distribution: expected {:?}, got {:?}",
            expected_chunks, fof_chunks
        );

        Ok(())
    }

    /*
        // The function "merge_bloom_filters" is no longer used
        #[test]
        fn test_merge_bloom_filters() {
            use roaring::RoaringBitmap;

            // initialize bf1: represents 2 datasets, abundance_number = 3, partitioned_bf_size = 2
            let mut bf1 = RoaringBitmap::new();  //110001,000100
            bf1.insert(0); // 1st dataset, 1st abundance level
            bf1.insert(1); // 1st dataset, 2nd abundance level
            bf1.insert(5); // 2nd dataset, 3rd abundance level
            bf1.insert(9); // 2nd dataset, 1rd abundance level, 2nd column

            // initialize bf2: represents 1 dataset, abundance_number = 3
            let mut bf2 = RoaringBitmap::new(); // 111,001
            bf2.insert(0); // 1st dataset, 1st abundance level
            bf2.insert(1); // 1st dataset, 2nd abundance level
            bf2.insert(2); // 1st dataset, 3rd abundance level
            bf2.insert(5); // 1st dataset, 3rd abundance level, 2nd column

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
            // second column: 000100 + 001 -> 000100001
            expected_bf.insert(12);
            expected_bf.insert();


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
    #[apply(findere_fixture)]
    fn test_build_and_query_index_with_chunks_and_merge(
        #[case] z: usize,
        random_directory: AutoRemoveDirectory,
    ) {
        let test_dir = random_directory.filename().to_str().unwrap();
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fof_path = format!("{}/fof.txt", test_dir);
        let filepaths = (1..=8)
            .map(|i| format!("{test_dir}/file{i}Q.fa"))
            .collect_vec();

        // write the fof
        {
            let mut fof_file = File::create(&fof_path).expect("Failed to create fof.txt");
            filepaths.iter().for_each(|path| {
                writeln!(fof_file, "{path}").expect("Failed to write to fof.txt");
            });
        }

        let data = [
            ("seq1", 30, "AAAAAAAAATAAAAAAAAAAAACACAGATCA"),
            ("seq2", 12, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT"),
            ("seq3", 1500, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA"),
            ("seq4", 1000, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG"),
            ("seq5", 450, "AAAAAAAAAAAAAAAAAAAAACACCCCTGGG"),
            ("seq6", 4, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA"),
            ("seq7", 45110, "AAAAAAAAAAAAAAAAAAAAACACCCCTGGG"),
            ("seq8", 75, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA"),
        ];

        for (file_path, (seq_id, ka_value, sequence)) in filepaths.iter().zip(data) {
            let mut file = File::create(file_path).expect("Failed to create FASTA file");
            writeln!(file, ">{seq_id} ka:f:{ka_value}").expect("Failed to write header");
            writeln!(file, "{sequence}").expect("Failed to write sequence");
        }

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 2,
            nb_color: 8,
            abundance_number: NonZero::new(255).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: false,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let query_results_path =
            create_build_query(parameters, chunks_size, threshold, filepaths, 0, test_dir);

        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [(">seq1 ka:f:30".to_string(), String::from("file1Q"), 29)];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
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
    #[apply(findere_fixture)]
    fn test_build_and_query_index_with_chunks_and_merge2(
        #[case] z: usize,
        random_directory: AutoRemoveDirectory,
    ) {
        let test_dir = random_directory.filename().to_str().unwrap();
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
            (
                &file3_path,
                ("seq3", 1500, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA"),
            ),
            (
                &file4_path,
                ("seq4", 1000, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG"),
            ),
            (
                &file5_path,
                ("seq5", 450, "AAAAAAAAAAAAAAAAAAAAACACCCCTGGG"),
            ),
            (&file6_path, ("seq6", 4, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA")),
        ] {
            let mut file = File::create(file_path).expect("Failed to create FASTA file");
            writeln!(file, ">{} ka:f:{}", seq_id, ka_value).expect("Failed to write header");
            writeln!(file, "{}", sequence).expect("Failed to write sequence");
        }

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024,
            partition_number: 2,
            nb_color: 6,
            abundance_number: NonZero::new(255).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: false,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let file_paths = vec![
            file1_path.clone(),
            file2_path.clone(),
            file3_path.clone(),
            file4_path.clone(),
            file5_path.clone(),
            file6_path.clone(),
        ];
        let query_results_path =
            create_build_query(parameters, chunks_size, threshold, file_paths, 5, test_dir);

        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [
            (">seq6 ka:f:4".to_string(), String::from("file3Q"), 1476),
            (">seq6 ka:f:4".to_string(), String::from("file6Q"), 4),
        ];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
    }

    //#[ignore]
    #[apply(findere_fixture)]
    fn test_build_and_query_index_with_chunks_and_merge3(
        #[case] z: usize,
        random_directory: AutoRemoveDirectory,
    ) {
        let test_dir = random_directory.filename().to_str().unwrap();
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
            (
                &file3_path,
                ("seq3", 1500, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA"),
            ),
            (
                &file4_path,
                ("seq4", 1000, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG"),
            ),
            (
                &file5_path,
                ("seq5", 60_000, "TAAAAAAAAAAAAAAAAAAAACACCCCTGGG"),
            ),
            (&file6_path, ("seq6", 4, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA")),
            (
                &file7_path,
                ("seq7", 1000, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG"),
            ),
            (
                &file8_path,
                ("seq8", 450, "AAAAAAAAAAAAAAAAAAAAACACCCCTGGG"),
            ),
            (&file9_path, ("seq9", 4, "GGGGAAAAAAAAAAAAAAAAAACAAAAAGAA")),
        ] {
            let mut file = File::create(file_path).expect("Failed to create FASTA file");
            writeln!(file, ">{} ka:f:{}", seq_id, ka_value).expect("Failed to write header");
            writeln!(file, "{}", sequence).expect("Failed to write sequence");
        }

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 2,
            nb_color: 9,
            abundance_number: NonZero::new(255).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: false,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let file_paths = vec![
            file1_path.clone(),
            file2_path.clone(),
            file3_path.clone(),
            file4_path.clone(),
            file5_path.clone(),
            file6_path.clone(),
            file7_path.clone(),
            file8_path.clone(),
            file9_path.clone(),
        ];
        let query_results_path =
            create_build_query(parameters, chunks_size, threshold, file_paths, 4, test_dir);

        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [(
            ">seq5 ka:f:60000".to_string(),
            String::from("file5Q"),
            59149,
        )];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
    }

    //#[ignore]
    #[apply(findere_fixture)]
    fn test_build_and_query_index_with_chunks_and_merge4(
        #[case] z: usize,
        random_directory: AutoRemoveDirectory,
    ) {
        let test_dir = random_directory.filename().to_str().unwrap();
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

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 256 * 256,
            partition_number: 2,
            nb_color: 5,
            abundance_number: NonZero::new(255).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: false,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let file_paths = vec![
            file1_path.clone(),
            file2_path.clone(),
            file3_path.clone(),
            file4_path.clone(),
            file5_path.clone(),
        ];
        let query_results_path =
            create_build_query(parameters, chunks_size, threshold, file_paths, 4, test_dir);

        let mut results = load_query_result_csv(query_results_path);
        let mut expected_results = [
            (
                ">seq5 ka:f:60000".to_string(),
                String::from("file5Q"),
                59149,
            ),
            (">seq6 ka:f:4".to_string(), String::from("file3Q"), 1476),
            (">seq6 ka:f:4".to_string(), String::from("file5Q"), 4),
        ];

        results.sort();
        expected_results.sort();

        assert_eq!(results, expected_results);
    }

    //#[ignore]
    #[apply(findere_fixture)]
    fn test_build_and_query_index_with_chunks_and_merge_longseq(
        #[case] z: usize,
        random_directory: AutoRemoveDirectory,
    ) {
        let test_dir = random_directory.filename().to_str().unwrap();
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
            (&file1_path, ("seq1", 30, "AAAAAAAAAAAAAAAAAAAAAACACAGATTT")),
            (&file2_path, ("seq2", 12, "AAAAAAAAAAAAAAAAAAAAACACAGATCAT")),
            (
                &file3_path,
                ("seq3", 1500, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA"),
            ),
            (
                &file4_path,
                ("seq4", 1000, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG"),
            ),
            (
                &file5_path,
                (
                    "seq5",
                    60_000,
                    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
                ),
            ),
            (&file6_path, ("seq6", 4, "AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA")),
            (
                &file7_path,
                ("seq7", 1000, "AAAAAAAAAAAAAAAAAAAAAACACCCCTGG"),
            ),
            (
                &file8_path,
                ("seq8", 450, "AAAAAAAAAAAAAAAAAAAAACACCCCTGGG"),
            ),
            (&file9_path, ("seq9", 4, "GGGGAAAAAAAAAAAAAAAAAACAAAAAGAA")),
        ] {
            let mut file = File::create(file_path).expect("Failed to create FASTA file");
            writeln!(file, ">{} ka:f:{}", seq_id, ka_value).expect("Failed to write header");
            writeln!(file, "{}", sequence).expect("Failed to write sequence");
        }

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 32,
            partition_number: 2,
            nb_color: 9,
            abundance_number: NonZero::new(255).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: false,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let file_paths = vec![
            file1_path.clone(),
            file2_path.clone(),
            file3_path.clone(),
            file4_path.clone(),
            file5_path.clone(),
            file6_path.clone(),
            file7_path.clone(),
            file8_path.clone(),
            file9_path.clone(),
        ];
        let query_results_path =
            create_build_query(parameters, chunks_size, threshold, file_paths, 1, test_dir);

        let results = load_query_result_csv(query_results_path);
        let expected_results = [
            //("seq5".to_string(), 4, 57549),
            (">seq2 ka:f:12".to_string(), String::from("file2Q"), 12),
        ];

        for (expected, actual) in expected_results.iter().zip(results.iter()) {
            assert_eq!(
                expected, actual,
                "Mismatch: expected {:?}, got {:?}",
                expected, actual
            );
        }
    }

    #[apply(findere_fixture)]
    fn test_color_graph(#[case] z: usize, random_directory: AutoRemoveDirectory) {
        let test_dir = random_directory.filename().to_str().unwrap();
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let fasta_path = format!("{}/test.fasta", test_dir);
        let output_path = format!("{}/colored_graph_output.fasta", test_dir);

        {
            let mut file = File::create(&fasta_path).expect("Failed to create test.fasta");
            writeln!(file, ">seq1 ka:f:29").expect("Failed to write header");
            writeln!(file, "TAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACAGATCA")
                .expect("Failed to write sequence");
            writeln!(file, ">seq2 ka:f:250").expect("Failed to write header");
            writeln!(file, "TTTTTAATGATCGATTTTTTTTTTTACCCCTGG").expect("Failed to write sequence");
        }

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 1,
            abundance_number: NonZero::new(255).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: false,
            canonical: false,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let _batch_size = 2;

        let mut index = Reindeer2::new(parameters, String::from(test_dir));

        let (_file_paths, index_dir) = index
            .build(vec![fasta_path.clone()], chunks_size, threshold, false)
            .expect("Failed to build index");

        let index_from_disk =
            Reindeer2::load_from_disk(&index_dir).expect("Failed to load index infos from disk");
        index_from_disk
            .query(
                &fasta_path,
                test_dir,
                &output_path,
                OutputFormat::Colored { normalized: None },
                0.5,
            )
            .expect("Failed to color graph");

        let expected_output = vec![
            ">seq1 ka:f:29 col:0:29",
            "TAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACAGATCA",
            ">seq2 ka:f:250 col:0:249",
            "TTTTTAATGATCGATTTTTTTTTTTACCCCTGG",
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
    }

    /// Asserts if the two contents are "equal", in the sense that:
    /// - their headers match
    /// - their **sorted** content match
    fn assert_equal_sorted_content_with_equal_header(content_a: &str, content_b: &str) {
        let mut lines_a = content_a.lines();
        let mut lines_b = content_b.lines();

        assert_eq!(lines_a.next().unwrap(), lines_b.next().unwrap());

        assert_eq!(
            lines_a.sorted().collect_vec(),
            lines_b.sorted().collect_vec()
        );
    }

    #[apply(findere_fixture)]
    fn test_output_rd1(#[case] z: usize, random_directory: AutoRemoveDirectory) {
        let test_dir = random_directory.filename().to_str().unwrap();
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let file1_path = String::from("tests/unit_tests_data/random_seq_with_revcomp.fa");
        let file2_path = String::from("tests/unit_tests_data/random_seq.fa");
        let file_paths = vec![file1_path, file2_path];

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 2,
            abundance_number: NonZero::new(256).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: false,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let mut index = Reindeer2::new(parameters, String::from(test_dir));
        let (_, index_dir) = index
            .build(file_paths.clone(), chunks_size, threshold, false)
            .expect("Failed to build index");

        let query_results_path = format!("{}/query_results.csv", index_dir);

        let index_from_disk = Reindeer2::load_from_disk(&index_dir).unwrap();
        index_from_disk
            .query(
                &file_paths[1],
                &index_dir,
                &query_results_path,
                OutputFormat::AbundanceMatrix {
                    format: MatrixFormat::Raw(None),
                },
                0.5,
            )
            .expect("Failed to query sequences");

        // Validate the results written to the query result file
        let actual = fs::read_to_string(&query_results_path).unwrap();
        let actual = actual.trim();

        let expected = String::from(
            "query\trandom_seq_with_revcomp\trandom_seq
header_0\t0-69:*\t0-69:1
header_1\t0-69:*\t0-69:1
header_2\t0-69:*\t0-69:1
header_3\t0-69:*\t0-69:1
header_4\t0-69:*\t0-69:1
shared_revcomp_with_other_test_file\t0-19:3\t0-19:10",
        );

        assert_equal_sorted_content_with_equal_header(&expected, actual);
    }

    #[apply(findere_fixture)]
    fn test_output_duplication(#[case] z: usize, random_directory: AutoRemoveDirectory) {
        let test_dir = random_directory.filename().to_str().unwrap();
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let file1_path = String::from("tests/unit_tests_data/random_seq_with_revcomp.fa");
        let file2_path = String::from("tests/unit_tests_data/duplication.fa");
        let file_paths = vec![file1_path, file2_path];

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 2,
            abundance_number: NonZero::new(256).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: false,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let mut index = Reindeer2::new(parameters, String::from(test_dir));
        let (_, index_dir) = index
            .build(file_paths.clone(), chunks_size, threshold, false)
            .expect("Failed to build index");

        let query_results_path = format!("{}/query_results.csv", index_dir);

        let index_from_disk = Reindeer2::load_from_disk(&index_dir).unwrap();
        index_from_disk
            .query(
                &file_paths[1],
                &index_dir,
                &query_results_path,
                OutputFormat::AbundanceMatrix {
                    format: MatrixFormat::Raw(None),
                },
                0.5,
            )
            .expect("Failed to query sequences");

        // Validate the results written to the query result file
        let actual = fs::read_to_string(&query_results_path).unwrap();
        let actual = actual.trim();

        let expected = String::from(
            "query\trandom_seq_with_revcomp\tduplication
header_0\t0-69:*\t0-69:1
header_0\t0-69:*\t0-69:1
header_0\t0-69:*\t0-69:1",
        );

        assert_equal_sorted_content_with_equal_header(&expected, actual);
    }

    #[apply(findere_fixture)]
    fn test_output_rd1_non_canonical(#[case] z: usize, random_directory: AutoRemoveDirectory) {
        let test_dir = random_directory.filename().to_str().unwrap();
        fs::create_dir_all(test_dir).expect("Failed to create test directory");

        let file1_path = String::from("tests/unit_tests_data/random_seq_with_revcomp.fa");
        let file2_path = String::from("tests/unit_tests_data/random_seq.fa");
        let file_paths = vec![file1_path, file2_path];

        let parameters = Parameters {
            k: 31,
            m: 15,
            bf_size: 1024 * 1024,
            partition_number: 4,
            nb_color: 2,
            abundance_number: NonZero::new(256).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: false,
            canonical: false,
            sampling_strategy: None,
            findere_z: z,
            capacity: wrong_capacity(),
        };
        let chunks_size = 128;
        let threshold = parameters.nb_color;

        let mut index = Reindeer2::new(parameters, String::from(test_dir));
        let (_, index_dir) = index
            .build(file_paths.clone(), chunks_size, threshold, false)
            .expect("Failed to build index");

        let query_results_path = format!("{}/query_results.csv", index_dir);

        let index_from_disk = Reindeer2::load_from_disk(&index_dir).unwrap();
        index_from_disk
            .query(
                &file_paths[1],
                &index_dir,
                &query_results_path,
                OutputFormat::AbundanceMatrix {
                    format: MatrixFormat::Raw(None),
                },
                0.5,
            )
            .expect("Failed to query sequences");

        // Validate the results written to the query result file
        let actual = fs::read_to_string(&query_results_path).unwrap();
        let actual = actual.trim();
        let expected = String::from(
            "query\trandom_seq_with_revcomp\trandom_seq
header_0\t0-69:*\t0-69:1
header_1\t0-69:*\t0-69:1
header_2\t0-69:*\t0-69:1
header_3\t0-69:*\t0-69:1
header_4\t0-69:*\t0-69:1
shared_revcomp_with_other_test_file\t0-19:*\t0-19:10",
        );

        assert_equal_sorted_content_with_equal_header(&expected, actual);
    }

    // TODO duplicate
    const MAX_SIZE_ROARING: u64 = 2u64.pow(32);
    /// Returns the minimum number of partitions required to index `nb_files` files with `nb_abundance_level` abundance levels and filters of size `bf_size`.
    #[must_use]
    pub fn get_number_of_partitions(
        nb_files: usize,
        nb_abundance_level: usize,
        bf_size: u64,
    ) -> usize {
        usize::try_from(
            (nb_abundance_level as u64 * nb_files as u64 * bf_size).div_ceil(MAX_SIZE_ROARING),
        )
        .expect("error in conversion")
    }

    #[apply(findere_fixture)]
    fn test_merge_multiple_indexes(
        #[case] z: usize,
        random_directory: AutoRemoveDirectory,
    ) -> io::Result<()> {
        let test_dir = random_directory.filename().to_str().unwrap();
        let index1_files_dir = format!("{}/index1_files", test_dir);
        let index2_files_dir = format!("{}/index2_files", test_dir);
        let index1_index_dir = format!("{}/index1_index", test_dir);
        let index2_index_dir = format!("{}/index2_index", test_dir);
        let merged_index_dir = format!("{}/merged_index", test_dir);
        let indexes_fof = format!("{}/indexes.txt", test_dir);

        fs::create_dir_all(&index1_files_dir)?;
        fs::create_dir_all(&index2_files_dir)?;
        fs::create_dir_all(&index1_index_dir)?;
        fs::create_dir_all(&index2_index_dir)?;

        let file1_path = format!("{}/file1.fa", index1_files_dir);
        {
            let mut file = File::create(&file1_path)?;
            writeln!(file, ">seq1 ka:f:30")?;
            writeln!(file, "TAGAAGGCGTGAGTCTTAGCT")?;
        }

        let file2_path = format!("{}/file2.fa", index2_files_dir);
        let file3_path = format!("{}/file3.fa", index2_files_dir);
        {
            let mut file = File::create(&file2_path)?;
            writeln!(file, ">seq2 ka:f:30")?;
            writeln!(file, "TAGAAGGCGTGAGTCTTAGCT")?;
        }
        {
            let mut file = File::create(&file3_path)?;
            writeln!(file, ">seq3 ka:f:30")?;
            writeln!(file, "TAGAAGGCGTGAGTCTTAGCT")?;
        }

        let capacity = 40;
        let bf_size = 1024;

        let partition_number = get_number_of_partitions(capacity, 255, bf_size);

        let mut parameters = Parameters {
            k: 21,
            m: 5,
            bf_size,
            partition_number,
            nb_color: 1,
            abundance_number: NonZero::new(255).unwrap(),
            abundance_min: 0,
            abundance_max: NonZero::new(65535).unwrap(),
            dense_option: false,
            canonical: true,
            sampling_strategy: None,
            findere_z: z,
            capacity,
        };

        let chunks_size = 128;
        let tolerated_number_of_zeros = 0;

        // for index1, color count = 1
        let index1_file_paths = vec![file1_path.clone()];
        let mut index = Reindeer2::new(parameters.clone(), String::from(&index1_index_dir));
        index.build(
            index1_file_paths,
            chunks_size,
            tolerated_number_of_zeros,
            false,
        )?;

        //index2, color count = 2
        parameters.nb_color = 2;
        let index2_file_paths = vec![file2_path.clone(), file3_path.clone()];
        let mut index = Reindeer2::new(parameters.clone(), String::from(&index2_index_dir));
        index.build(
            index2_file_paths,
            chunks_size,
            tolerated_number_of_zeros,
            false,
        )?;

        {
            let mut fof = File::create(&indexes_fof)?;
            writeln!(fof, "{}", index1_index_dir)?;
            writeln!(fof, "{}", index2_index_dir)?;
        }

        merge_multiple_indexes(&indexes_fof, &merged_index_dir)
            .expect("Failed to merge the test indexes");

        let merged_index = Reindeer2::load_from_disk(&merged_index_dir)
            .expect("Failed to read the merged index metadata)");
        assert_eq!(
            merged_index.parameters.nb_color, 3,
            "Expected merged color count to be 3, got {}",
            merged_index.parameters.nb_color
        );

        assert_eq!(
            merged_index.indexed_file_names,
            vec!["file1", "file2", "file3"],
        );

        // check that each partition file in the merged index exists
        let total_capacity: u64 = (0..NB_FILE_IN_AN_INDEX)
            .map(|i| {
                if i < parameters.partition_number {
                    let path = format!(
                        "{}/partition_bloom_filters_group{}.bin",
                        merged_index.index_dir, i
                    );
                    let file = File::open(path).unwrap();
                    let mut reader = BufReader::new(file);
                    let metadata = tar_get::get_metadata(&mut reader).unwrap();
                    metadata.get_nb_objects()
                } else {
                    0
                }
            })
            .sum();
        assert_eq!(total_capacity, parameters.partition_number as u64);

        Ok(())
    }
}
