mod cli;
mod reindeer2;

use clap::Parser;
use rand::Rng;
use std::io::{self};
use std::time::Instant;

use crate::reindeer2::{
    build_index, build_index_muset, explore_muset_dir, query_index, read_fof_file,
};
use cli::Cli;

fn main() -> io::Result<()> {
    let args = Cli::parse();

    let max_threads = args.threads;
    let debug = args.debug;

    match args.command {
        cli::Command::Index(args) => {
            let input = args.input;
            let kmer = args.kmer;
            let minimizer = args.minimizer;
            let partitions = args.partitions;
            let bloomfilter = args.bloomfilter;
            let abundance = args.abundance;
            let abundance_max = args.abundance_max;
            let dense_option = args.dense;
            let output_dir = args.output_dir.unwrap_or_else(|| {
                format!("PACAS_index_{}", rand::thread_rng().gen::<u64>()) // Generate a unique directory name
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
                    panic!("ERROR : With the '--dense' option set to 'true', the k-mer size must be <= 32.")
                }
                if abundance > 255 {
                    println!("WARNING : the abundance granularity exceeds the requirements of the '--dense' (<256). The abundance granularity is now set to 255.");
                }
            } else {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(max_threads)
                    .build_global()
                    .unwrap();
            }
            let minimizer = if kmer < minimizer {
                println!("WARNING : the minimizer size '{}' exceeds the k-mer size '{}'. The minimiser size is now set to '{}'", minimizer, kmer, kmer);
                kmer
            } else {
                minimizer
            };
            println!();

            let tolerated_number_of_zeros = 0;

            let start_time = Instant::now();

            if false {
                //muset_option {

                let (unitigs_file, matrix_file, color_nb) = explore_muset_dir(&input);

                build_index_muset(
                    unitigs_file,
                    matrix_file,
                    kmer,
                    minimizer,
                    bf_size,
                    partitions,
                    color_nb,
                    abundance,
                    abundance_max,
                    &output_dir,
                    dense_option,
                    tolerated_number_of_zeros,
                    debug,
                )?;
            } else {
                // read the file of files  and extract file paths and color count
                let (file_paths, color_nb) = read_fof_file(&input)?;

                // run the index construction process: build and fill BFs per partitions and in chunks, serialize, merge chunks
                build_index(
                    file_paths,
                    kmer,
                    minimizer,
                    bf_size,
                    partitions,
                    color_nb,
                    abundance,
                    abundance_max,
                    &output_dir,
                    dense_option,
                    tolerated_number_of_zeros,
                    debug,
                )?;
            }

            println!("Indexing complete in {:.2?}", start_time.elapsed());
        }

        cli::Command::Query(args) => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(max_threads)
                .build_global()
                .unwrap();

            let fasta_file = args.fasta;
            let index_dir = args.index;
            let color_graph = args.color;
            let coverage = args.coverage;
            let normalize_option = args.normalize;

            println!("Index directory: {}", index_dir);

            let mut query_output = format!("{}/query_results.csv", index_dir);
            if color_graph {
                query_output = format!("{}/colored_graph.fa", index_dir);
            }

            let start_time = Instant::now();
            query_index(
                &fasta_file,
                &index_dir,
                &query_output,
                color_graph,
                normalize_option,
                coverage,
            )
            .expect("Failed to query sequences");

            println!("Query complete in {:.2?}", start_time.elapsed());
        } // "merge" => {
          //     // argument= path to a fof + output file
          //     // let indexes_fof = matches
          //     //     .get_one::<String>("indexes")
          //     //     .expect("Required argument: indexes (file-of-index directories)");
          //     // let output_dir = matches.get_one::<String>("output");
          //     // merge_multiple_indexes(indexes_fof, output_dir.as_deref().map(|x| x.as_str()))?;
          // }
          // _ => {
          //     eprintln!("Invalid mode: {}. Use 'index' or 'query'.", mode);
          //     // eprintln!("Invalid mode: {}. Use 'index', 'query', or 'merge'.", mode);
          // }
    }

    Ok(())
}
