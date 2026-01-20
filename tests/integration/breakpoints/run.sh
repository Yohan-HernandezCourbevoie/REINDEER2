#!/bin/bash

source tests/integration/breakpoints/check_file_contents.sh

cargo build --quiet

reindeer2="./target/debug/reindeer2"

LOCAL_FOLDER="tests/integration/breakpoints"

INPUT_FOF="$LOCAL_FOLDER/variations.fof"

QUERY_INPUT="$LOCAL_FOLDER/query.fa"
QUERY_OUTPUT="$LOCAL_FOLDER/results.csv"
EXPECTED_QUERY_OUTPUT="$LOCAL_FOLDER/expected.csv"

RUST_LOG=warn $reindeer2 index --input $INPUT_FOF -k 31 -o integration_test_index
RUST_LOG=warn $reindeer2 query --fasta $QUERY_INPUT --index ./integration_test_index --output-format abundance-matrix-raw --breakpoints 0.3 --output $QUERY_OUTPUT

check_file_contents $QUERY_OUTPUT $EXPECTED_QUERY_OUTPUT
is_same=$?
echo $is_same

rm -r ./integration_test_index 
rm $QUERY_OUTPUT

# Compare
if [ $is_same -eq 0 ]; then
    exit 0
else
    echo "Test failed"
    exit 1
fi


