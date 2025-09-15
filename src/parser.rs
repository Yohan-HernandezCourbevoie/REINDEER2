use clap::{Arg, ArgAction, ArgMatches, Command};

// TODO discuss utiliser un struct ici
pub fn parse_args() -> ArgMatches {
    Command::new("REINDEER 2")
        .bin_name("Reindeer2")
        .version("1.0")
        .author("Camille Marchet, Yohan Hernandez--Courbevoie")
        .about("REINDEER2\nk-mer quantification in large collections of samples")
        .after_help("Example:\n  $ Reindeer2 --mode index --input test_files/fof.txt --kmer 31 --output-dir ../index_test\n  $ Reindeer2 --mode query --fasta test_files/file1Q.fa --index ../index_test")
        .disable_version_flag(true)
        .disable_help_flag(true)
        .arg(
            Arg::new("mode")
                .short('M')
                .long("mode")
                .value_name("MODE")
                .help("Mode: 'index' to build the index or 'query' to query a built index\n")
                .required(true),
        )
        .arg(
            Arg::new("input")
                .short('I')
                .long("input")
                .value_name("INPUT")
                .help("For 'index' mode (required): By default, a file of files where each line is the path to a FASTA file (Logan format). With (-u, --muset) set, a path to a muset
                output directory."),
        )
        .arg(
            Arg::new("kmer")
                .short('k')
                .long("kmer")
                .value_name("SIZE")
                .help("For 'index' mode (required): Sets the k-mer size"),
        )
        .arg(
            Arg::new("minimizer")
                .short('m')
                .long("minimizer")
                .value_name("MINSIZE")
                .help("For 'index' mode: Sets the minimizer size (default: 15)"),
        )
        .arg(
            Arg::new("partitions")
                .short('p')
                .long("partitions")
                .value_name("MINSIZE")
                .help("For 'index' mode: Sets the number of partitions (default: 512)"),
        )
        .arg(
            Arg::new("bloomfilter")
                .short('b')
                .long("bloomfilter")
                .value_name("BF")
                .help("For 'index' mode: Sets the Bloom filter size in log2 scale (default: 32)"),
        )
        .arg(
            Arg::new("abundance")
                .short('a')
                .long("abundance")
                .value_name("ABUND")
                .help("For 'index' mode: Sets the abundance granularity (default: 255)"),
        )
        .arg(
            Arg::new("abundance_max")
                .short('A')
                .long("abundance-max")
                .value_name("ABUND_MAX")
                .help("For 'index' mode: Sets the maximal abundance to take into account (default: 65024)"),
        )
        .arg(
            Arg::new("dense")
                .short('d')
                .long("dense")
                .value_name("DENSE")
                .action(ArgAction::Set)
                .help("For 'index' mode: If set, allows to index dense k-mers - i.e. shared k-mers among datasets - more efficiently, at the cost of higher RAM consumption, \
                limited parameters (k-mer size <= 32, number of abundance levels <= 255) and limited multithreading (default: false)")
        )
        // .arg(
        //     Arg::new("threshold")
        //         .short('T')
        //         .long("threshold")
        //         .value_name("THRESHOLD")
        //         .action(ArgAction::Set)
        //         .help("For 'index' mode (with --dense true): Minimum proportion of occurrence in samples for a k-mer to be considered omnipresent, 0 < T <= 1 (default: 1.0)")
        // )
        // .arg(
        //     Arg::new("muset")
        //         .short('u')
        //         .long("muset")
        //         .value_name("MUSET")
        //         .action(ArgAction::Set)
        //         .help("For 'index' mode: If set, the index takes as input the output directory of Muset, containing at least 'unitigs.fa' and 'unitigs.abundance.mat' (default: false)")
        // )
        .arg(
            Arg::new("output_dir")
                .short('o')
                .long("output-dir")
                .value_name("OUT")
                .help("For 'index' mode: Sets the index output directory\n"),
        )
        .arg(
            Arg::new("fasta")
                .short('f')
                .long("fasta")
                .value_name("FILE")
                .help("For 'query' mode (required): Path to the FASTA file containing query sequences"),
        )
        .arg(
            Arg::new("index")
                .short('i')
                .long("index")
                .value_name("DIR")
                .help("For 'query' mode (required): Path to the directory containing the prebuilt index"),
        )
        .arg(
            Arg::new("color")
                .short('c')
                .long("color")
                .value_name("COLOR")
                .help("For 'query' mode: bool for coloring a graph instead of regular query (default: false)"),
        )
        .arg(
            Arg::new("normalize")
                .short('n')
                .long("normalize")
                .value_name("NORMALIZE")
                .help("For 'query' mode: bool for normalizing abundances based on sequencing depth estimates (default: false)"),
        )
        .arg(
            Arg::new("coverage")
                .short('C')
                .long("coverage-min")
                .value_name("COVERAGE")
                .action(ArgAction::Set)
                .help("For 'query' mode : Minimum proportion of kmers that must be present in the query sequence in order to propose an abundance value, 0 < C <= 1 (default: 0.5)")
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .value_name("THREADS")
                .action(ArgAction::Set)
                .help("Define a number of threads available (default: 1)")
        )
        .arg(
            Arg::new("debug")
                .long("debug")
                .value_name("DEBUG")
                .help("Show additionnal informations for debugging purposes\n"),
        )
        .arg(
            Arg::new("version")
                .short('V')
                .long("version")
                .global(true)
                .action(ArgAction::Version)
                .help("Print version"),
        )
        .arg(
            Arg::new("help")
                .short('h')
                .long("help")
                .global(true)
                .action(ArgAction::Help)
                .help("Print help"),
        )
        .get_matches()
}
