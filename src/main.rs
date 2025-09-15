mod parser;
mod reindeer2;

use crate::parser::parse_args;
use crate::reindeer2::{
    build_index, build_index_muset, explore_muset_dir, query_index, read_fof_file,
};

use std::io::{self};
use std::time::Instant;

use rand::Rng;

fn main() -> io::Result<()> {
    let matches = parse_args();

    let mode = matches
        .get_one::<String>("mode")
        .expect("Required argument for mode (index or query)");
    let max_threads = matches
        .get_one::<String>("threads")
        .map(|s| s.parse::<usize>().expect("Invalid number of threads"))
        .unwrap_or(1); // default number of threads

    match mode.as_str() {
        "index" => {
            // PARAMETERS
            let input = matches
                .get_one::<String>("input")
                .expect("Required argument for INPUT in index mode");

            let kmer = matches
                .get_one::<String>("kmer")
                .expect("Required argument for k-mer size")
                .parse::<usize>()
                .expect("Invalid k-mer size");

            let minimizer = matches
                .get_one::<String>("minimizer")
                .map(|s| s.parse::<usize>().expect("Invalid minimizer size"))
                .unwrap_or(15); // default size

            let partitions = matches
                .get_one::<String>("partitions")
                .map(|s| s.parse::<usize>().expect("Invalid number of partitions"))
                .unwrap_or(512); // default number of partitions

            let bloomfilter = matches
                .get_one::<String>("bloomfilter")
                .map(|s| s.parse::<usize>().expect("Invalid Bloom filter size"))
                .unwrap_or(32); // default size

            let bf_size = 1u64 << bloomfilter; // Bloom filter size as a power of 2

            let abundance = matches
                .get_one::<String>("abundance")
                .map(|s| s.parse::<usize>().expect("Invalid abundance number"))
                .unwrap_or(255); // default abundance levels

            let abundance_max = matches
                .get_one::<String>("abundance_max")
                .map(|s| s.parse::<u16>().expect("Invalid maximal abundance"))
                .unwrap_or(65024);

            let dense_option = matches
                .get_one::<String>("dense")
                .map(|s| s.parse::<bool>().expect("Invalid color option"))
                .unwrap_or(false);

            // let muset_option = matches
            //     .get_one::<String>("muset")
            //     .map(|s| s.parse::<bool>().expect("Invalid color option"))
            //     .unwrap_or(false);

            // TODO add threads

            let output_dir = match matches.get_one::<String>("output_dir") {
                Some(v) => v,
                None => {
                    // maybe use a UUID ?
                    let mut rng = rand::thread_rng();
                    let dir_seed: u64 = rng.gen();
                    &format!("PACAS_index_{}", dir_seed) // Generate a unique directory name
                }
            };

            let debug = matches
                .get_one::<String>("debug")
                .map(|s| s.parse::<bool>().expect("Invalid debug option"))
                .unwrap_or(false);

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

                let (unitigs_file, matrix_file, color_nb) = explore_muset_dir(input);

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
                    output_dir,
                    dense_option,
                    tolerated_number_of_zeros,
                    debug,
                )?;
            } else {
                // read the file of files  and extract file paths and color count
                let (file_paths, color_nb) = read_fof_file(input)?;

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
                    output_dir,
                    dense_option,
                    tolerated_number_of_zeros,
                    debug,
                )?;
            }

            println!("Indexing complete in {:.2?}", start_time.elapsed());
        }

        "query" => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(max_threads)
                .build_global()
                .unwrap();

            let fasta_file = matches
                .get_one::<String>("fasta")
                .expect("Required argument for FASTA file in query mode");

            let index_dir = matches
                .get_one::<String>("index")
                .expect("Required argument for index directory in query mode");

            let color_graph = matches
                .get_one::<String>("color")
                .map(|s| s.parse::<bool>().expect("Invalid color option"))
                .unwrap_or(false);

            let coverage = matches
                .get_one::<String>("coverage")
                .map(|s| s.parse::<f32>().expect("Invalid coverage threshold value"))
                .unwrap_or(0.5); // default size

            let normalize_option = matches
                .get_one::<String>("normalize")
                .map(|s| s.parse::<bool>().expect("Invalid normalize option"))
                .unwrap_or(false);

            println!("Index directory: {}", index_dir);

            let mut query_output = format!("{}/query_results.csv", index_dir);
            if color_graph {
                query_output = format!("{}/colored_graph.fa", index_dir);
            }

            let start_time = Instant::now();
            query_index(
                fasta_file,
                index_dir,
                &query_output,
                color_graph,
                normalize_option,
                coverage,
            )
            .expect("Failed to query sequences");

            println!("Query complete in {:.2?}", start_time.elapsed());
        }
        // "merge" => {
        //     // argument= path to a fof + output file
        //     // let indexes_fof = matches
        //     //     .get_one::<String>("indexes")
        //     //     .expect("Required argument: indexes (file-of-index directories)");
        //     // let output_dir = matches.get_one::<String>("output");
        //     // merge_multiple_indexes(indexes_fof, output_dir.as_deref().map(|x| x.as_str()))?;
        // }
        _ => {
            eprintln!("Invalid mode: {}. Use 'index' or 'query'.", mode);
            // eprintln!("Invalid mode: {}. Use 'index', 'query', or 'merge'.", mode);
        }
    }

    Ok(())
}
