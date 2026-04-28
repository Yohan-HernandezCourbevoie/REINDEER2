#!/bin/bash

##### variable definitions #####
reindeer2="./target/debug/reindeer2"  # we must define `reindeer2` as it is not installed yet 
LOCAL_FOLDER="tests/system_testing/colored"
INPUT_FOF="$LOCAL_FOLDER/datasets/colored.fof"
QUERY_INPUT="$LOCAL_FOLDER/query.fa"
QUERY_OUTPUT="$LOCAL_FOLDER/results.fa"
EXPECTED_QUERY_OUTPUT="$LOCAL_FOLDER/expected.fa"

##### use of REINDEER2 #####
cargo build --quiet  # build REINDEER2
# we can configure the log level, setting to warn only prints the warnings
RUST_LOG=warn $reindeer2 index --input $INPUT_FOF -k 31 -o integration_test_index --no-sort-files-by-size
RUST_LOG=warn $reindeer2 query --fasta $QUERY_INPUT --index ./integration_test_index --output-format colored --output $QUERY_OUTPUT

##### ensure the example works #####
python3 "$LOCAL_FOLDER/files_equal.py" $QUERY_OUTPUT $EXPECTED_QUERY_OUTPUT
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


