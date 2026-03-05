mod cli;
mod memory_measure;
mod overflow_detection;

use clap::Parser;
use rand::Rng;
use std::io::{self};
use std::num::NonZero;
use std::time::Instant;

use cli::Cli;
use cli::OutputFormatCli;
use overflow_detection::{get_min_number_of_files, get_number_of_partitions};
use reindeer2::reindeer2::{
    merge_multiple_indexes, read_fof_file, BreakpointsNormalize, MatrixFormat, OutputFormat,
    Parameters, Reindeer2, SamplingStrategy,
};

use crate::cli::{IndexArgs, InfosArgs, MergeArgs, QueryArgs};

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
fn main() -> io::Result<()> {
    let args = Cli::parse();
    let env = env_logger::Env::default().default_filter_or("trace");
    env_logger::init_from_env(env);
    let max_threads = args.threads;

    match args.command {
        cli::Command::Index(IndexArgs {
            input,
            kmer,
            minimizer,
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
        }) => {
            let dense_option = dense;
            let canonical = !stranded;
            let output_dir = output_dir.unwrap_or_else(|| {
                format!("RD2_index_{}", rand::rng().random::<u64>()) // Generate a unique directory name
            });
            // let muset_option = args.muset;
            // TODO add threads

            let bf_size = 1u64 << bloomfilter; // Bloom filter size as a power of 2

            // CHECKS
            if dense_option {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(1)
                    .build_global()
                    .unwrap();
                if kmer > 32 {
                    panic!(
                        "ERROR : With the '--dense' option set to 'true', the k-mer size must be <= 32."
                    )
                }
            } else {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(max_threads)
                    .build_global()
                    .unwrap();
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

            let tolerated_number_of_zeros = 0;

            let start_time = Instant::now();

            // if false {
            //     //muset_option {

            //     let (unitigs_file, matrix_file, color_nb) = explore_muset_dir(&input);

            //     build_index_muset(
            //         unitigs_file,
            //         matrix_file,
            //         kmer,
            //         minimizer,
            //         bf_size,
            //         partitions,
            //         color_nb,
            //         abundance,
            //         abundance_max,
            //         &output_dir,
            //         dense_option,
            //         tolerated_number_of_zeros,
            //         canonical,
            //         debug,
            //     )?;
            // } else {
            // read the file of files  and extract file paths and color count
            let (file_paths, color_nb) = read_fof_file(&input)?;

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
            };
            let mut index = Reindeer2::new(parameters, output_dir);
            index.build(file_paths, chunks_size, tolerated_number_of_zeros)?;
            // }

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

        cli::Command::Query(QueryArgs {
            fasta,
            index,
            output_format,
            normalize,
            output,
            coverage,
            breakpoints,
        }) => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(max_threads)
                .build_global()
                .unwrap();

            let fasta_file = fasta;
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
            let index = Reindeer2::load_from_disk(&index_dir)
                .expect("should have been able to load index infos from disk");
            index
                .query(
                    &fasta_file,
                    &index_dir,
                    &query_output,
                    output_format,
                    coverage,
                )
                .expect("Failed to query sequences");
            log::info!("Results written to {}", query_output);
            log::info!("Query complete in {:.2?}", start_time.elapsed());
        }

        cli::Command::Merge(MergeArgs {
            file_of_indexes,
            output_dir,
        }) => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(max_threads)
                .build_global()
                .unwrap();

            let output_dir = output_dir.unwrap_or_else(|| {
                format!("RD2_index_{}", rand::rng().random::<u64>()) // Generate a unique directory name
            });

            let start_time = Instant::now();
            merge_multiple_indexes(&file_of_indexes, &output_dir)
                .expect("Failed to merge the given indexes.");

            log::info!("Query complete in {:.2?}", start_time.elapsed());
        }
        cli::Command::Infos(InfosArgs { index }) => {
            let index_dir = index;
            let index = Reindeer2::load_from_disk(&index_dir)
                .expect("should have been able to load index infos from disk");
            let infos = index.get_index_infos();

            println!("{infos}");
        }
    }

    Ok(())
}
