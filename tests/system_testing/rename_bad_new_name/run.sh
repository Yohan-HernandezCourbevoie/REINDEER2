#!/bin/bash

##### variable definitions #####
reindeer2=$(realpath "./target/debug/reindeer2")  # we must define `reindeer2` as it is not installed yet
LOCAL_FOLDER=$(dirname $0)
INPUT_FOF="datasets/fof.fof"
QUERY_INPUT="query.fa"
QUERY_OUTPUT="results.tsv"
EXPECTED_QUERY_OUTPUT="expected.tsv"
EXPECTED_ERROR="datasetB already appears in the list of indexed files"

cd $LOCAL_FOLDER

cargo build --quiet
RUST_LOG=warn $reindeer2 index --input $INPUT_FOF -k 31 -o integration_test_index
RD2OUTPUT=$(RUST_LOG=warn $reindeer2 rename --index integration_test_index --old-name datasetA --new-name datasetB 2>&1)

rm -r ./integration_test_index 

# Compare
if [[ "$RD2OUTPUT" == *"$EXPECTED_ERROR"* ]]; then
    exit 0
else
    echo "Test failed"
    exit 1
fi