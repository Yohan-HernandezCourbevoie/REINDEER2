#!/bin/bash

##### variable definitions #####
reindeer2=$(realpath "./target/debug/reindeer2")  # we must define `reindeer2` as it is not installed yet
LOCAL_FOLDER=$(dirname $0)
INPUT_FOF_0="datasets/matrix_raw_0.fof"
INPUT_FOF_1="datasets/matrix_raw_1.fof"
INPUT_FOF_INDEXES="datasets/matrix_raw_index.fof"
QUERY_INPUT="query.fa"
QUERY_OUTPUT="results.tsv"
EXPECTED_QUERY_OUTPUT="expected.tsv"

cd $LOCAL_FOLDER

cargo build --quiet
RUST_LOG=warn $reindeer2 index --threads 7 --input $INPUT_FOF_0 -k 31 --nb-file-capacity 2 -o integration_test_index_0
RUST_LOG=warn $reindeer2 index --threads 7 --input $INPUT_FOF_1 -k 31 --nb-file-capacity 2 -o integration_test_index_1
RUST_LOG=warn $reindeer2 merge --threads 7 --file-of-indexes $INPUT_FOF_INDEXES -o integration_test_index
RUST_LOG=warn $reindeer2 query --fasta $QUERY_INPUT --index ./integration_test_index --output-format matrix-raw --output $QUERY_OUTPUT -C 0

python3 files_equal.py $QUERY_OUTPUT $EXPECTED_QUERY_OUTPUT
is_same=$?

rm -r ./integration_test_index_0 
rm -r ./integration_test_index_1 
rm -r ./integration_test_index 
rm $QUERY_OUTPUT

# Compare
if [ $is_same -eq 0 ]; then
    exit 0
else
    echo "Test failed"
    exit 1
fi