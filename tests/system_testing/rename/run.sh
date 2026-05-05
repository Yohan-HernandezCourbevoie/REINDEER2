#!/bin/bash

##### variable definitions #####
reindeer2=$(realpath "./target/debug/reindeer2")  # we must define `reindeer2` as it is not installed yet
LOCAL_FOLDER=$(dirname $0)
INPUT_FOF="datasets/fof.fof"
QUERY_INPUT="query.fa"
QUERY_OUTPUT="results.tsv"
EXPECTED_QUERY_OUTPUT="expected.tsv"
EXPECTED_OUTPUT="datasetA was renamed to datasetC"

cd $LOCAL_FOLDER

cargo build --quiet
RUST_LOG=warn $reindeer2 index --input $INPUT_FOF -k 31 -o integration_test_index
RUST_LOG=warn $reindeer2 infos integration_test_index > output_before_rename.txt
RD2OUTPUT=$(RUST_LOG=warn $reindeer2 rename --index integration_test_index --old-name datasetA --new-name datasetC)
RUST_LOG=warn $reindeer2 infos integration_test_index > output_after_rename.txt
RUST_LOG=warn $reindeer2 query --fasta $QUERY_INPUT --index ./integration_test_index --output-format matrix-raw --output $QUERY_OUTPUT

python3 check_before.py output_before_rename.txt
is_same_before_rename=$?
python3 check_after.py output_after_rename.txt
is_same_after_rename=$?
python3 files_equal.py $QUERY_OUTPUT $EXPECTED_QUERY_OUTPUT
is_same_query_output=$?

rm -r ./integration_test_index 
rm -r output_before_rename.txt
rm -r output_after_rename.txt

# Compare
if [[ "$RD2OUTPUT" == *"$EXPECTED_OUTPUT"* ]]; then
    exit 0
else
    echo "Test failed"
    exit 1
fi

if [ $is_same_before_rename -eq 0 ]; then
    exit 0
else
    echo "Test failed"
    exit 1
fi

if [ $is_same_after_rename -eq 0 ]; then
    exit 0
else
    echo "Test failed"
    exit 1
fi

if [ $is_same_query_output -eq 0 ]; then
    exit 0
else
    echo "Test failed"
    exit 1
fi