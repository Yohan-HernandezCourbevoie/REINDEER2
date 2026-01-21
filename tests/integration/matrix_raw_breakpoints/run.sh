#!/bin/bash

reindeer2="./target/debug/reindeer2"
LOCAL_FOLDER="tests/integration/matrix_raw_breakpoints"
INPUT_FOF="$LOCAL_FOLDER/datasets/matrix_raw_breakpoints.fof"
QUERY_INPUT="$LOCAL_FOLDER/query.fa"
QUERY_OUTPUT="$LOCAL_FOLDER/results.tsv"
EXPECTED_QUERY_OUTPUT="$LOCAL_FOLDER/expected.tsv"

cargo build --quiet

RUST_LOG=warn $reindeer2 index --input $INPUT_FOF -k 31 -o integration_test_index
RUST_LOG=warn $reindeer2 query --fasta $QUERY_INPUT --index ./integration_test_index --output-format abundance-matrix-raw --breakpoints 0.3 --output $QUERY_OUTPUT

python3 "$LOCAL_FOLDER/check_files_are_similar.py" $QUERY_OUTPUT $EXPECTED_QUERY_OUTPUT
is_same=$?

rm -r ./integration_test_index 
rm $QUERY_OUTPUT

# Compare
if [ $is_same -eq 0 ]; then
    exit 0
else
    echo "Test failed"
    exit 1
fi


