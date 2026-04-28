#!/bin/bash

reindeer2="./target/debug/reindeer2"
LOCAL_FOLDER="tests/system_testing/matrix_raw_sampled"
INPUT_FOF="$LOCAL_FOLDER/datasets/matrix_raw.fof"
QUERY_INPUT="$LOCAL_FOLDER/query.fa"
QUERY_OUTPUT="$LOCAL_FOLDER/results.tsv"
EXPECTED_QUERY_OUTPUT="$LOCAL_FOLDER/expected.tsv"

cargo build --quiet

EXPECTED_QUERY_OUTPUT_0="$LOCAL_FOLDER/expected_0.tsv"
EXPECTED_QUERY_OUTPUT_1="$LOCAL_FOLDER/expected_1.tsv"
EXPECTED_QUERY_OUTPUT_2="$LOCAL_FOLDER/expected_2.tsv"
EXPECTED_QUERY_OUTPUT_3="$LOCAL_FOLDER/expected_3.tsv"
EXPECTED_QUERY_OUTPUT_4="$LOCAL_FOLDER/expected_4.tsv"
EXPECTED_QUERY_OUTPUT_5="$LOCAL_FOLDER/expected_5.tsv"

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

