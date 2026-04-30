#!/bin/sh

# Get repository root
REPO_ROOT=$(git rev-parse --show-toplevel)

# Run test script at repo root
TEST_DIR="$REPO_ROOT/tests/system_testing"
cd "$REPO_ROOT"

pids=""
failed=""

for file in "$TEST_DIR"/*/run.sh; do
    [ -f "$file" ] || continue
    tmpfile=$(mktemp)
    "$file" >"$tmpfile" 2>&1 &
    pids="$pids $! $file $tmpfile"
done

# Wait for all jobs and collect failures
set -- $pids
while [ $# -ge 3 ]; do
    pid=$1 file=$2 tmpfile=$3
    shift 3
    if ! wait "$pid"; then
        failed="$failed $file"
        printf '\n=== FAILED: %s ===\n' "$file"
        cat "$tmpfile"
    fi
    rm -f "$tmpfile"
done

if [ -n "$failed" ]; then
    printf '\nFailed scripts:%s\n' "$(printf '\n  %s' $failed)"
    exit 1
fi

exit 0