#!/bin/sh

# Get repository root
REPO_ROOT=$(git rev-parse --show-toplevel)

# Run test script at repo root
TEST_DIR="$REPO_ROOT/tests/system_testing"
cd "$REPO_ROOT"
for file in "$TEST_DIR"/*/run.sh; do
    [ -f "$file" ] || continue   # skip if no match
    "$file" || { printf '%s\n' "$file"; exit 1; }            # stop if any script fails
done

exit 0