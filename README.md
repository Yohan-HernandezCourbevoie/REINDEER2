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

In the following examples, the tool's command will use its installation name `Reindeer2` but it can be replaced by `cargo run --` when used from the REINDEER2 folder.


## Usage

The generale use of REINDEER 2 is divided in to steps : index building and abundance query.

### Index

For the **index** mode, the mandatory parameters are the file of files (a plain text file where each line represented a unitigs file) and the size of the k-mers to be indexed.

`Reindeer2 --mode index --input file_of_files.txt --kmer 31`


General parameters:
- `-o, --output-dir` an output directory for the index
- `-a, --abundance` the abundance granularity (number of levels or discretized abundance values)
- `-A, --abundance-max` the maximal abundance to take into account
- `-d, --dense` (true/false) allows to index dense k-mers - shared k-mers among datasets - more efficiently (default: false)
<!-- - `-u, --muset` (true/false) the index takes as input the output directory of Muset, containing at least 'unitigs.fa' and 'unitigs.abundance.mat' (default: false) -->
- `-t, --threads` the maximal number of threads used (default: 1)

Advanced parameters: 
- `-b, --bloomfilter` the Bloom filter size in log2 scale
- `-m, --minimizer` the minimizer size
- `-p, --partitions` the number of partitions



### Query

For **query** mode, the parameters are the FASTA file containing the sequence(s) to be queried and the index directory.

`Reindeer2 --mode query --fasta sequences_query.fa --index ~/index_directory`

Optional parameters:
- `-c, --color` (true/false) is used to annotate the input file with abundances rather than producing the standard output file (as showed in the examples below) (default: false)
- `-n, --normalize` (true/false) parameter allows to normalize abundances based on sequencing depth estimates. The calculation is _normalized\_abundance = raw\_abundance / number\_of\_kmers\_in\_the\_dataset * 1\_000\_000_ (default: false)
- `-C, --coverage-min` minimum proportion of kmers that must be present in the query sequence in order to propose an abundance value

#### CSV file with --color false (default)

This option outputs a CSV file (header included) with the following structure : `<Sequence_header>,<Color>,<Median_abundance>`

#### FASTA file with --color true

This option outputs the FASTA file given in query annotated with the abundances of all indexed files.

## Example

To illustrate how the tool works, a small example is available. The commands are launched from the REINDEER2 main directory.

#### INDEX
How to build the index:
```
Reindeer2 --mode index --input test_files/fof.txt --kmer 31 --output-dir ../index_test
```

#### QUERY (results: CSV)
With the command:
```
Reindeer2 --mode query --fasta test_files/file1Q.fa --index ../index_test
cat ../index_test/query_results.csv
```
is expected the result:
```
header,file,abundance
>seq1 ka:f:30,0,29
>seq2 ka:f:30,0,29
>seq3 ka:f:1500,0,1470
>seq3 ka:f:1500,1,4
```


#### QUERY (results: colored graph FASTA)
With the command:
```
Reindeer2 --mode query --fasta test_files/file1Q.fa --index ../index_test --color true 
cat ../index_test/colored_graph.fa
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

