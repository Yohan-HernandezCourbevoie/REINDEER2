use clap::{Args, Parser, Subcommand, ValueEnum};

#[derive(Parser, Debug)]
// #[command(disable_version_flag = true, disable_help_flag = true)]// TODO discuss
#[command(
    after_help = "Example:\n  $ Reindeer2 --mode index --input test_files/fof.txt --kmer 31 --output-dir ../index_test\n  $ Reindeer2 --mode query --fasta test_files/file1Q.fa --index ../index_test",
    about, long_about = None
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,

    /// Define a number of threads available (default: 1)
    #[arg(
        short = 't',
        long = "threads",
        value_name = "THREADS",
        default_value_t = 1
    )]
    pub threads: usize,

    /// Show additional information for debugging purposes
    #[arg(long = "debug")]
    pub debug: bool,
    // /// Print version
    // #[arg(short = 'V', long = "version", global = true, action = ArgAction::Version)]
    // pub version: bool,

    // /// Print help
    // #[arg(short = 'h', long = "help", global = true, action = ArgAction::Help)]
    // pub help: bool,
}

// TODO name the modes
#[derive(Subcommand, Debug)]
pub enum Command {
    /// Index file
    #[clap(alias = "b")]
    Index(IndexArgs),
    /// Query file
    #[clap(alias = "d")]
    Query(QueryArgs),
}

#[derive(Args, Debug)]
pub struct IndexArgs {
    /// By default, a file of files where each line is the path to a FASTA file (Logan format).
    /// With (-u, --muset) set, a path to a muset output directory.
    #[arg(short = 'I', long = "input", value_name = "INPUT")]
    pub input: String,

    /// Sets the k-mer size
    #[arg(short = 'k', long = "kmer", value_name = "SIZE")]
    pub kmer: usize,

    /// Sets the minimizer size (default: 15)
    #[arg(
        short = 'm',
        long = "minimizer",
        value_name = "MINSIZE",
        default_value_t = 15
    )]
    pub minimizer: usize,

    // TODO check default is printed any ?
    /// Sets the number of partitions (default: 512)
    #[arg(
        short = 'p',
        long = "partitions",
        value_name = "MINSIZE",
        default_value_t = 512
    )]
    pub partitions: usize,

    /// Sets the Bloom filter size in log2 scale (default: 32)
    #[arg(
        short = 'b',
        long = "bloomfilter",
        value_name = "BF",
        default_value_t = 32
    )]
    pub bloomfilter: usize,

    /// Sets the abundance granularity (default: 255)
    #[arg(
        short = 'a',
        long = "abundance",
        value_name = "ABUND",
        default_value_t = 255
    )]
    pub abundance: usize,

    /// Sets the maximal abundance to take into account (default: 65024)
    #[arg(
        short = 'A',
        long = "abundance-max",
        value_name = "ABUND_MAX",
        default_value_t = 65024
    )]
    pub abundance_max: u16,

    /// If set, allows to index dense k-mers - i.e. shared k-mers among datasets - more efficiently,
    /// at the cost of higher RAM consumption, limited parameters (k-mer size <= 32, number of abundance levels <= 255)
    /// and limited multithreading (default: false)
    #[arg(short = 'd', long = "dense", default_value_t = false)]
    pub dense: bool,

    /// Sets the index output directory (default: random name in the form of PACAS_index_)
    #[arg(short = 'o', long = "output-dir", value_name = "OUT")]
    pub output_dir: Option<String>,

    /// Use non-canonical version of k-mers (default: false)
    #[arg(long, default_value_t = false)]
    pub stranded: bool,
}

#[derive(Copy, Clone, Debug, PartialEq, ValueEnum)]
pub enum OutputFormat {
    Colored,
    NormalizedMedian,
    Median,
    Raw,
}

#[derive(Args, Debug)]
pub struct QueryArgs {
    /// Path to the FASTA file containing query sequences
    #[arg(short = 'f', long = "fasta", value_name = "FILE")]
    pub fasta: String,

    /// Path to the directory containing the prebuilt index
    #[arg(short = 'i', long = "index", value_name = "DIR")]
    pub index: String,

    /// Precision of the format output
    #[arg(long, value_enum, default_value_t = OutputFormat::Median)]
    pub output_format: OutputFormat,

    /// Precision of the format output
    #[arg(short, long)]
    pub output: Option<String>,

    // TODO discuss that 0.5 default value
    /// Minimum proportion of kmers that must be present in the query sequence in order to propose an abundance value, 0 < C <= 1 (default: 0.5)
    #[arg(
        short = 'C',
        long = "coverage-min",
        value_name = "COVERAGE",
        default_value_t = 0.5
    )]
    pub coverage: f32,
}
