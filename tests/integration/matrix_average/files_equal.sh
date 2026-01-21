#!/bin/bash

files_equal() {
    local file1="$1"
    local file2="$2"
    local line_num=0
    local diff_lines=()

    # Check if files exist and are readable
    if [[ ! -f "$file1" || ! -f "$file2" ]]; then
        echo "Error: One or both files do not exist." >&2
        return 2
    fi

    # Compare files line by line
    while IFS= read -r line1 <&3 && IFS= read -r line2 <&4; do
        ((line_num++))
        if [[ "$line1" != "$line2" ]]; then
            diff_lines+=("$line_num: $line1")
            diff_lines+=("$line_num: $line2")
            break
        fi
    done 3<"$file1" 4<"$file2"

    # Check if one file has more lines than the other
    if [[ -n "$line1" || -n "$line2" ]]; then
        if [[ -z "$line1" ]]; then
            diff_lines+=("EOF: $file1")
            diff_lines+=("$line_num: $line2")
        else
            diff_lines+=("$line_num: $line1")
            diff_lines+=("EOF: $file2")
        fi
    fi

    # Log differences if any
    if [[ ${#diff_lines[@]} -gt 0 ]]; then
        echo "Files differ:"
        for line in "${diff_lines[@]}"; do
            echo "$line"
        done
        return 1
    else
        echo "Files are equal."
        return 0
    fi
}



