use clap::{Args, Parser, Subcommand, ValueEnum};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,

    /// Define a number of threads available
    #[arg(short, long, value_name = "THREADS", default_value_t = 1)]
    pub threads: usize,
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
    /// Merge option
    #[clap(alias = "m")]
    Merge(MergeArgs),
}

#[derive(Args, Debug)]
pub struct IndexArgs {
    /// A file of files where each line is the path to a multi-FASTA file of unitigs (Logan format).
    #[arg(short = 'I', long = "input", value_name = "INPUT")]
    pub input: String,

    /// Sets the k-mer size
    #[arg(short = 'k', long = "kmer", value_name = "SIZE")]
    pub kmer: usize,

    /// Sets the minimizer size
    #[arg(short, long, value_name = "MINSIZE", default_value_t = 15)]
    pub minimizer: usize,

    /// Sets the number of partitions
    #[arg(short, long, value_name = "PARTS", default_value_t = 512)]
    pub partitions: usize,

    /// Sets the Bloom filter size in log2 scale
    #[arg(short, long, value_name = "BF", default_value_t = 32)]
    pub bloomfilter: usize,

    /// Sets the abundance granularity
    #[arg(short, long, value_name = "ABUND", default_value_t = 255)]
    pub abundance: usize,

    /// Sets the maximal abundance to take into account
    #[arg(short = 'A', long, value_name = "ABUND_MAX", default_value_t = 65024)]
    pub abundance_max: u16,

    /// Sets the number of datasets treated at a time, affecting RAM consumption
    #[arg(short, long, value_name = "chunksize", default_value_t = 128)]
    pub chunks_size: usize,

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
pub enum OutputFormatCli {
    Colored,
    Median,
    AbundanceMatrixRaw,
    AbundanceMatrixMedian,
    AbundanceMatrixAverage,
}

#[derive(Args, Debug)]
pub struct QueryArgs {
    /// Path to the FASTA file containing query sequences
    #[arg(short, long, value_name = "FILE")]
    pub fasta: String,

    /// Path to the directory containing the prebuilt index
    #[arg(short, long, value_name = "DIR")]
    pub index: String,

    /// Precision of the format output
    #[arg(long, value_enum, default_value_t = OutputFormatCli::Median)]
    pub output_format: OutputFormatCli,

    /// Normalize output (default: false)
    #[arg(long, num_args(0..=1), value_parser = clap::value_parser!(u64), default_missing_value = "1000000")]
    pub normalize: Option<u64>,

    /// Path to the output file
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

    /// Allows to detect the (non ordered) set of breakpoints in the output.
    #[arg(short, long)]
    pub breakpoints: Option<f64>,
}

#[derive(Args, Debug)]
pub struct MergeArgs {
    /// A file of indexes where each line is the path an index directory.
    #[arg(short = 'f', long = "file-of-indexes", value_name = "FILE_OF_INDEXES")]
    pub file_of_indexes: String,

    /// Sets the index output directory (default: random name in the form of PACAS_index_)
    #[arg(short = 'o', long = "output-dir", value_name = "OUT")]
    pub output_dir: Option<String>,
}
