#!/bin/bash

##### variable definitions #####
reindeer2=$(realpath "./target/debug/reindeer2")  # we must define `reindeer2` as it is not installed yet
LOCAL_FOLDER=$(dirname $0)
INPUT_FOF="datasets/fof.fof"
QUERY_INPUT="query.fa"
QUERY_OUTPUT="results.fa"
EXPECTED_QUERY_OUTPUT="expected.fa"

cd $LOCAL_FOLDER

##### use of REINDEER2 #####
cargo build --quiet  # build REINDEER2
# we can configure the log level, setting to warn only prints the warnings
RUST_LOG=warn $reindeer2 index --input $INPUT_FOF -k 31 -o integration_test_index --no-sort-files-by-size
RUST_LOG=warn $reindeer2 query --fasta $QUERY_INPUT --index ./integration_test_index --output-format colored --output $QUERY_OUTPUT

##### ensure the example works #####
python3 files_equal.py $QUERY_OUTPUT $EXPECTED_QUERY_OUTPUT
is_same=$?

##### cleanup #####
rm -r ./integration_test_index 
rm $QUERY_OUTPUT

##### Report error if necessary #####
if [ $is_same -eq 0 ]; then
    exit 0
else
    echo "Test failed"
    exit 1
fi


