mod cli;
mod memory_measure;
mod overflow_detection;

use clap::Parser;
use log::warn;
use rand::Rng;
use std::io::{self};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::time::Instant;

use cli::Cli;
use cli::OutputFormatCli;
use overflow_detection::{get_min_number_of_files, get_number_of_partitions};
use reindeer2::{
    BreakpointsNormalize, BuildArgs, MatrixFormat, OutputFormat, Parameters, Reindeer2,
    SamplingStrategy, merge::merge_multiple_indexes, read_fof_file,
};

use crate::cli::{IndexArgs, InfosArgs, MergeArgs, QueryArgs, RenameArgs, ResumeIndexationArgs};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

impl OutputFormatCli {
    fn to_output_format(self, normalized: Option<u64>, breakpoints: Option<f64>) -> OutputFormat {
        match (self, normalized, breakpoints) {
            // raw matrix
            (OutputFormatCli::MatrixRaw, Some(normalized), None) => OutputFormat::AbundanceMatrix {
                format: MatrixFormat::Raw(Some(BreakpointsNormalize::Normalize(normalized))),
            },
            (OutputFormatCli::MatrixRaw, None, Some(penalty)) => OutputFormat::AbundanceMatrix {
                format: MatrixFormat::Raw(Some(BreakpointsNormalize::Breakpoints(penalty))),
            },
            (OutputFormatCli::MatrixRaw, None, None) => OutputFormat::AbundanceMatrix {
                format: MatrixFormat::Raw(None),
            },
            (OutputFormatCli::MatrixRaw, Some(_), Some(_)) => {
                panic!("Cannot compute both breakpoints and normalized abundance.")
            }

            // matrix average
            (OutputFormatCli::MatrixAverage, normalized, None) => OutputFormat::AbundanceMatrix {
                format: MatrixFormat::Average { normalized },
            },
            (OutputFormatCli::MatrixAverage, _, Some(_)) => {
                panic!("Cannot compute breakpoints from a matrix with mode average.")
            }

            // matrix median
            (OutputFormatCli::MatrixMedian, normalized, None) => OutputFormat::AbundanceMatrix {
                format: MatrixFormat::Median { normalized },
            },
            (OutputFormatCli::MatrixMedian, _normalized, Some(_)) => {
                panic!("Cannot compute breakpoints from a matrix with mode median.")
            }

            // colored
            (OutputFormatCli::Colored, normalized, None) => OutputFormat::Colored { normalized },
            (OutputFormatCli::Colored, _, Some(_)) => {
                panic!("Cannot compute breakpoints from colored graph.")
            }

            // median
            (OutputFormatCli::Median, normalized, None) => OutputFormat::Median { normalized },
            (OutputFormatCli::Median, _, Some(_)) => {
                panic!("Cannot compute breakpoints from median output.")
            }
        }
    }
}

fn validate_sampling_strategy(
    kmer_sampling: Option<u64>,
    minimizer_sampling: Option<u64>,
) -> Option<SamplingStrategy> {
    match (kmer_sampling, minimizer_sampling) {
        (Some(_), Some(_)) => {
            panic!("cannot compute sampling from both kmer sampling and minimizer sampling");
        }
        (Some(kmer_sampling_factor), None) => Some(SamplingStrategy::KmerSampling {
            last_bits_to_zero: kmer_sampling_factor,
        }),
        (None, Some(minimizer_sampling_factor)) => Some(SamplingStrategy::MinimizerSampling {
            last_bits_to_zero: minimizer_sampling_factor,
        }),
        (None, None) => None,
    }
}

#[cfg(feature = "self-destruct")]
fn validate_fail(
    chunk_explode_at_step: Option<usize>,
    merge_explode_at_step: Option<usize>,
) -> Option<reindeer2::FailIndexation> {
    match (chunk_explode_at_step, merge_explode_at_step) {
        (None, None) => None,
        (Some(a), None) => Some(reindeer2::FailIndexation::Chunk(a)),
        (None, Some(b)) => Some(reindeer2::FailIndexation::Merge(b)),
        (Some(_), Some(_)) => {
            panic!("cannot error twice");
        }
    }
}

#[cfg(feature = "self-destruct")]
fn self_destruct_warn() {
    log::warn!(
        "You are using the \"self-destruct\" feature. Only use it for testing. Consider recompiling the tool without this feature."
    );
    println!(
        "Warning: you are using the \"self-destruct\" feature. Only use it for testing. Consider recompiling the tool without this feature."
    );
}

fn main() -> io::Result<()> {
    let args = Cli::parse();
    let env = env_logger::Env::default().default_filter_or("trace");
    env_logger::init_from_env(env);

    match args.command {
        cli::Command::Index(IndexArgs {
            input,
            kmer,
            minimizer,
            threads,
            nb_file_capacity,
            bloomfilter,
            abundance,
            abundance_min,
            abundance_max,
            chunks_size,
            dense,
            stranded,
            output_dir,
            kmer_sampling,
            minimizer_sampling,
            allow_count_right_after_angle_bracket,
            findere,
            no_sort_files_by_size,
            #[cfg(feature = "self-destruct")]
            chunk_explode_at_step,
            #[cfg(feature = "self-destruct")]
            merge_explode_at_step,
        }) => {
            #[cfg(feature = "self-destruct")]
            self_destruct_warn();

            let sort_files_by_size = !no_sort_files_by_size;
            let dense_option = dense;
            let canonical = !stranded;
            let output_dir = output_dir.unwrap_or_else(|| {
                format!("RD2_index_{}", rand::rng().random::<u64>()) // Generate a unique directory name
            });

            let bf_size = 1u64 << bloomfilter; // Bloom filter size as a power of 2

            // CHECKS
            if dense_option {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(1)
                    .build_global()
                    .expect("should have been able to set up the threads (maybe the setup function was called twice ?)");
                if kmer > 32 {
                    panic!(
                        "ERROR : With the '--dense' option set to 'true', the k-mer size must be <= 32."
                    )
                }
            } else {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build_global()
                    .expect("should have been able to set up the threads (maybe the setup function was called twice ?)");
            }
            let abundance = if abundance > 255 && dense_option {
                log::warn!(
                    "WARNING : the abundance granularity exceeds the requirements of the '--dense' (<256). The abundance granularity is now set to 255."
                );
                255
            } else {
                abundance
            };
            let abundance =
                NonZero::new(abundance).expect("abundance granularity should be positive");

            let abundance_max =
                NonZero::new(abundance_max).expect("max abundance should be positive");

            assert!(
                abundance_min < abundance_max.get(),
                "the minimum abundance ({abundance_min}) should not be smaller than the max abundance ({abundance_max})"
            );

            let minimizer = if kmer < minimizer {
                log::warn!(
                    "WARNING : the minimizer size '{}' exceeds the k-mer size '{}'. The minimiser size is now set to '{}'",
                    minimizer,
                    kmer,
                    kmer
                );
                kmer
            } else {
                minimizer
            };

            let sampling_strategy = validate_sampling_strategy(kmer_sampling, minimizer_sampling);

            // test that findere is not activated at the same time as sampling
            let findere_z = match (findere, sampling_strategy) {
                (Some(_), Some(_)) => panic!("sampling and findere are currently not supported"),
                (Some(z), None) => z, // user asked for findere without sampling
                (None, Some(_)) => 0, // user asked sampling and not findere
                (None, None) => cli::DEFAULT_Z, // user did'nt asked anything, let's provide findere's default
            };

            assert!(
                kmer > findere_z,
                "Fatal error: findere's z parameter cannot be higher than k. We recommand (k-z) > 16."
            );

            if (kmer - findere_z) <= 16 {
                warn!(
                    "Indexing with k = {} and z = {}. A high false positive rate is expected. Using (k-z) > 16 is recommanded.",
                    kmer, findere_z
                );
                println!(
                    "Warning: with current chosen values (k = {}, findere's z = {}), the index might have a lot of false positives. We recommand using using (k-z) > 16.",
                    kmer, findere_z
                );
            }

            let tolerated_number_of_zeros = 0;

            let start_time = Instant::now();

            // read the file of files  and extract file paths and color count
            let (file_paths, color_nb) = read_fof_file(&input).unwrap_or_else(|err| {
                panic!("should have been able to read the input file {input} ({err})")
            });

            let nb_files = get_min_number_of_files(&file_paths, nb_file_capacity);
            let partitions = get_number_of_partitions(nb_files, abundance.get(), bf_size);

            // run the index construction process: build and fill BFs per partitions and in chunks, serialize, merge chunks
            let parameters = Parameters {
                bf_size,
                partition_number: partitions,
                k: kmer,
                m: minimizer,
                nb_color: color_nb,
                abundance_number: abundance,
                abundance_min,
                abundance_max,
                dense_option,
                canonical,
                sampling_strategy: validate_sampling_strategy(kmer_sampling, minimizer_sampling),
                findere_z,
                capacity: nb_files,
                #[cfg(feature = "self-destruct")]
                fail: validate_fail(chunk_explode_at_step, merge_explode_at_step),
            };
            let output_dir = PathBuf::from(output_dir);
            let mut index = Reindeer2::new(parameters, output_dir);
            let input_path = Path::new(&input);
            let input_file_name = input_path
                .file_name()
                .expect("impossible to extract the name of the input file")
                .to_str()
                .expect("the input file's name is not in UTF-8")
                .to_string();
            let build_args = BuildArgs {
                file_of_file_name: input_file_name,
                sort_files_by_size,
                chunks_size,
                threshold: tolerated_number_of_zeros,
                allow_count_right_after_angle_bracket,
            };
            index.build(build_args, file_paths)?;

            log::info!("Indexing complete in {:.2?}", start_time.elapsed());
            if cfg!(unix) {
                // linux returns memory in kb, no idea if all platforms do the same
                // just print on linux to prevent being inaccurate
                log::info!(
                    "Peak memory usage: {} kibibytes",
                    memory_measure::format_int_with_spaces(memory_measure::get_max_rss())
                );
            }
        }
        cli::Command::ResumeIndexation(ResumeIndexationArgs {
            partial_index_directory,
            #[cfg(feature = "self-destruct")]
            chunk_explode_at_step,
            #[cfg(feature = "self-destruct")]
            merge_explode_at_step,
        }) => {
            #[cfg(feature = "self-destruct")]
            self_destruct_warn();
            let partial_index_directory = PathBuf::from(&partial_index_directory);
            let build_args = Reindeer2::load_build_args(&partial_index_directory).unwrap_or_else(|| {
                log::error!("No incomplete index found. Please check the index you are trying to resume exists and is incomplete.");
        panic!(
"Should have been able to load the build arguments from disk. Please check the index you are trying to resume exists and is incomplete."        )
    });

            let mut index = Reindeer2::load_from_disk(partial_index_directory)
                .expect("should have been able to load index infos from disk");
            // index.restart_build(build_args)?;
        }

        cli::Command::Query(QueryArgs {
            fasta,
            index,
            output_format,
            threads,
            normalize,
            output,
            coverage,
            breakpoints,
        }) => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .expect("should have been able to set up the threads (maybe the setup function was called twice ?)");

            let fasta_file = PathBuf::from(fasta);
            let index_dir = index;
            let output_format = output_format.to_output_format(normalize, breakpoints);
            let query_output = match output {
                Some(output) => output,
                None => match output_format {
                    OutputFormat::Colored { normalized: _ } => {
                        format!("{}_colored_graph.fa", index_dir)
                    }
                    OutputFormat::Median { normalized: _ } => {
                        format!("{}_query_results.csv", index_dir)
                    }

                    _ => format!("{}_query_results.tsv", index_dir),
                },
            };

            log::info!("Index directory: {}", index_dir);

            let start_time = Instant::now();
            let index_dir = PathBuf::from(index_dir);
            let query_output = Path::new(&query_output);
            let index = Reindeer2::load_from_disk(index_dir.clone())
                .expect("should have been able to load index infos from disk");
            index
                .query(
                    &fasta_file,
                    &index_dir,
                    query_output,
                    output_format,
                    coverage,
                )
                .expect("Failed to query sequences");
            log::info!("Results written to {}", query_output.display());
            log::info!("Query complete in {:.2?}", start_time.elapsed());
        }

        cli::Command::Merge(MergeArgs {
            file_of_indexes,
            output_dir,
            threads,
        }) => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .expect("should have been able to set up the threads (maybe the setup function was called twice ?)");

            let output_dir = output_dir.unwrap_or_else(|| {
                format!("RD2_index_{}", rand::rng().random::<u64>()) // Generate a unique directory name
            });
            let output_dir = PathBuf::from(output_dir);
            let file_of_indexes = PathBuf::from(file_of_indexes);

            let start_time = Instant::now();
            merge_multiple_indexes(&file_of_indexes, &output_dir)
                .expect("Failed to merge the given indexes.");

            log::info!("Merge complete in {:.2?}", start_time.elapsed());
        }
        cli::Command::Infos(InfosArgs { index }) => {
            let index_dir = PathBuf::from(index);
            let index = Reindeer2::load_from_disk(index_dir)
                .expect("should have been able to load index infos from disk");
            let infos = index.get_index_infos();

            println!("{infos}");
        }
        cli::Command::Rename(RenameArgs {
            index,
            old_name,
            new_name,
        }) => {
            let index_dir = PathBuf::from(index);
            let mut index = Reindeer2::load_from_disk(index_dir)
                .expect("should have been able to load index infos from disk");
            let outcome = index
                .rename(&old_name, new_name.clone())
                .expect("should have been able to write updated index to disk");

            use reindeer2::ReplaceOutcome;
            match outcome {
                ReplaceOutcome::NotFound => {
                    log::error!("{old_name} was not found");
                }
                ReplaceOutcome::WouldAppearMoreThanOnce => {
                    log::error!("{new_name} already appears in the list of indexed files");
                }
                ReplaceOutcome::Replaced => {
                    println!("{old_name} was renamed to {new_name}");
                }
            };
        }
    }

    Ok(())
}
