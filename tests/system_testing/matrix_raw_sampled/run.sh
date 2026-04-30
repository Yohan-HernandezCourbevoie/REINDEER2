#!/bin/bash

##### variable definitions #####
reindeer2=$(realpath "./target/debug/reindeer2")  # we must define `reindeer2` as it is not installed yet
LOCAL_FOLDER=$(dirname $0)
INPUT_FOF="datasets/fof.fof"
QUERY_INPUT="query.fa"
QUERY_OUTPUT="results.tsv"
EXPECTED_QUERY_OUTPUT="expected.tsv"

EXPECTED_QUERY_OUTPUT_0="expected_0.tsv"
EXPECTED_QUERY_OUTPUT_1="expected_1.tsv"
EXPECTED_QUERY_OUTPUT_2="expected_2.tsv"
EXPECTED_QUERY_OUTPUT_3="expected_3.tsv"
EXPECTED_QUERY_OUTPUT_4="expected_4.tsv"
EXPECTED_QUERY_OUTPUT_5="expected_5.tsv"

cd $LOCAL_FOLDER

cargo build --quiet

set -euo pipefail

for sampling in 0 1 2 3 4; do
    expected_var="EXPECTED_QUERY_OUTPUT_${sampling}"
    expected="${!expected_var}"

    RUST_LOG=warn "$reindeer2" index --input "$INPUT_FOF" -k 31 --kmer-sampling "$sampling" -o integration_test_index --no-sort-files-by-size
    RUST_LOG=warn "$reindeer2" query --fasta "$QUERY_INPUT" --index ./integration_test_index --output-format matrix-raw --output "$QUERY_OUTPUT"

    python3 "$LOCAL_FOLDER/files_equal.py" "$QUERY_OUTPUT" "$expected"

    rm -r ./integration_test_index
    rm "$QUERY_OUTPUT"
done

