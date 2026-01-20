#!/bin/bash

cargo run index --input test_files/variations.fof -k 31 -o integration_test_index
cargo run query --fasta test_files/varations_query.fa --index ./integration_test_index --output-format abundance-matrix-raw --breakpoints 0.3 --output integration_test.csv
rm -r ./integration_test_index 
# cat integration_test.csv  # expect set 0,20,70



# Expected values
expected="0,20,70"

# Read actual values
actual=$(<integration_test.csv)
rm integration_test.csv

# Convert to sorted lists
sorted_expected=$(echo "$expected" | tr ',' '\n' | sort -n | tr '\n' ',')
sorted_actual=$(echo "$actual" | tr ',' '\n' | sort -n | tr '\n' ',')

# Remove trailing commas
sorted_expected=${sorted_expected%,}
sorted_actual=${sorted_actual%,}

# Compare
if [ "$sorted_actual" = "$sorted_expected" ]; then
    echo "Test passed"
    exit 1
else
    echo "Test failed"
    echo "Expected: $sorted_expected"
    echo "Got: $sorted_actual"
    exit 1
fi


