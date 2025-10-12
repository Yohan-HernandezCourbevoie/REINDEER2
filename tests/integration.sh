# TODO make a lib
# TODO once the lib is done, make real integration tests
cargo run index --input test_files/variations.fof -k 31 -o integration_test_index
cargo run query --fasta test_files/varations_query.fa --index ./integration_test_index --output-format abundance-matrix --breakpoints 0.3 --output integration_test.csv
rm -r ./integration_test_index 
echo ""
echo "=== results ==="
cat integration_test.csv  # expect set 0,20,70
rm integration_test.csv