# REINDEER 2

REINDEER 2 is an efficient and scalable k-mer abundance index.

## Pre-processing

As REINDEER 2 only indexes unitig files, a pre-processing step is necessary. This can be done using assembly tools such as [GGCAT](https://github.com/algbio/ggcat).

For public datasets, the sequencing files may already have been processed into unitigs by the [Logan project](https://github.com/IndexThePlanet/Logan). These files are freely available and can be downloaded by following [these steps](https://github.com/IndexThePlanet/Logan/blob/main/Accessions.md).   

## Installation

### Requirements

- cargo >= 1.81.0
- rustc >= 1.81.0

### Compilation

Clone the repository :

```
git clone https://github.com/Yohan-HernandezCourbevoie/REINDEER2.git 
```

Then build :

```
cd REINDEER2 && RUSTFLAGS="-C target-cpu=native" cargo build --release
```

Alternatively, the tool can be installed globally in the system using :

```
cd REINDEER2 && cargo install --path .
```

In the following examples, the tool's command uses its installation name `reindeer2`, but it can be replaced with `cargo run --` when used from the REINDEER2 folder.


## Usage

The general use of REINDEER 2 is divided into two steps: index building and abundance query.

### Index

For the **index** mode, the mandatory parameters are the file of files (a plain text file where each line represents a unitigs file) and the size of the k-mers to be indexed.

`reindeer2 index --input <file_of_files.txt> --kmer <k>`


General parameters:
- `-o, --output-dir` an output directory for the index (default: random name in the form of RD2_index_)
- `-a, --abundance` the abundance granularity (number of levels or discretized abundance values) (default: 255) (see ![the relation between abundance levels and approximation ratio](/doc/levels_and_approximation.pdf))
- `-A, --abundance-max` the maximal abundance to take into account
- `-d, --dense` (true/false) allows to index dense k-mers - shared k-mers among datasets - more efficiently (default: false)
<!-- - `-u, --muset` (true/false) the index takes as input the output directory of Muset, containing at least 'unitigs.fa' and 'unitigs.abundance.mat' (default: false) -->
- `-t, --threads` the maximal number of threads used (default: 1)
- `--stranded` use non-canonical version of k-mers (default: false)
- `-c, --chunks-size` number of datasets treated at a time, affecting RAM consumption (default: 128)
- `--abundance-min` minimum abundance for a k-mer to be indexed (default: 0, i.e., index all k-mers)

Advanced parameters: 
- `-b, --bloomfilter` the Bloom filter size in log2 scale (default: 32)
- `-m, --minimizer` the minimizer size (default: 15)
- `--nb-file-capacity` maximum number of files in the index. Default: reserve space only for the indexed files.
- `--no-sort-files-by-size` index the files without sorting them by their size (sorting files makes the indexation faster when multithreaded). Default: false (i.e., sort the input to speed up indexation). 

<!-- The value of the number of partitions is important when building large indexes. -->
\* More information on how to choose the right parameters [here](doc/parameters.md).


### Query
Warning: REINDEER 2's query results do not respect the order given in the indexed file of file. Please parse the query output rather than assuming the query will repect the order of the indexed file of file.


For **query** mode, the parameters are the FASTA file containing the sequence(s) to be queried and the index directory.

`reindeer2 query --fasta <sequences_query.fa> --index <index_directory>`

- `-f, --output-format` allows to change the output format. Supported formats are:
    - `median` (default): for each color, returns the median of k-mer abundance per read
    - `colored`: annotate the input file with abundances rather than producing the standard output file (as shown in the examples below)
    - `matrix-raw`: for each color, for each queried sequence, write a tsv containing the abundance of each k-mer (similar to REINDEER 1).
    - `matrix-median`: for each color, for each queried sequence, write a tsv containing the median of k-mers.
    - `matrix-average`: for each color, for each queried sequence, write a tsv containing the average of k-mers.
- `--breakpoints <penalty>`: Reindeer2 will apply the `PELT` algorithm to detect breakpoints in the abundances of k-mers. Reindeer will then report the position of such breakpoints in the query. This option is only available if the output format is `matrix-raw`. **Warning:** using this options significantly slows down the query.
- `-t, --threads` the maximal number of threads used (default: 1)
- `--normalize <N>`: normalize abundances based on sequencing depth estimates. The calculation is _normalized\_abundance = raw\_abundance / number\_of\_kmers\_in\_the\_dataset * N_. No normalization by default. If `--normalize` is passed without an argument, `N` defaults to 1\_000\_000. This option is incompatible with `--breakpoints`.
- `-C, --coverage-min` minimum proportion of kmers that must be present in the query sequence in order to propose an abundance value

### Infos

The **infos** command prints informations about an index.

`reindeer2 infos <index_path>`


<!--
#### CSV file with --color false (default)

This option outputs a CSV file (header included) with the following structure : `<Sequence_header>,<Color>,<Median_abundance>`

#### FASTA file with --color true

This option outputs the FASTA file given in query annotated with the abundances of all indexed files.
 -->

## Example

To illustrate how REINDEER2 works, examples are available in the folder `tests/system_testing`. See `tests/system_testing/colored/run.sh` for a commented example. 

The commands are launched from the REINDEER2 main directory.

#### INDEX
How to build the index:
```
reindeer2 index --input test_files/fof.txt --kmer 31 --output-dir ../index_test
```

#### QUERY (results: CSV)
With the command:
```
reindeer2 query --fasta test_files/file1Q.fa --index ../index_test
cat ../index_test_query_results.csv
```
is expected the result:
```
header,file,abundance
>seq1 ka:f:30,file1Q,29
>seq2 ka:f:30,file1Q,29
>seq3 ka:f:1500,file1Q,1470
>seq3 ka:f:1500,file2Q,4
```


#### QUERY (results: colored graph FASTA)
With the command:
```
reindeer2 query --fasta test_files/file1Q.fa --index ../index_test --output-format colored 
cat ../index_test_colored_graph.fa
```
is expected the result:
```
>seq1 ka:f:30 col:0:29 col:1:0
AAAAAAAAAAAAAAAAAAAAAACACAGATCA
>seq2 ka:f:30 col:0:29 col:1:0
AAAAAAAAAAAAAAAAAAAAACACAGATCAT
>seq3 ka:f:1500 col:0:1470 col:1:4
AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA
```

#### INFOS
With the command:
```
reindeer2 index --input test_files/fof.txt --kmer 31 --output-dir ../index_test
reindeer2 infos ../index_test
```
is expected the result:
```
index directory: "../index_test"
parameters: Parameters {
    bf_size: 4294967296,
    partition_number: 510,
    k: 31,
    m: 15,
    nb_color: 2,
    abundance_number: 255,
    abundance_min: 0,
    abundance_max: 65024,
    dense_option: false,
    canonical: true,
    sampling_strategy: None,
    findere_z: 4,
    capacity: 2,
}
indexed filenames and k-mers: [
  ("file1Q", 0),
  ("file2Q", 1000),
]
```

