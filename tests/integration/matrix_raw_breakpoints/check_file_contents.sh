#!/bin/bash

# Compare two files ignoring:
#  - first lines are equal
#  - first word of each remaining line is equal
#  - second word of each remaining line contains the same numbers
check_file_contents() {
    local file1="$1"
    local file2="$2"

    # Check if files exist and are readable
    if [[ ! -r "$file1" || ! -r "$file2" ]]; then
        echo "Error: One or both files are not readable."
        return 2
    fi

    # Read first line of each file
    local first_line1 first_line2
    first_line1=$(head -n 1 "$file1")
    first_line2=$(head -n 1 "$file2")

    # Check if first lines are equal
    if [[ "$first_line1" != "$first_line2" ]]; then
        return 1
    fi

    # Read remaining lines (skip first line)
    local lines1 lines2
    mapfile -t lines1 < <(tail -n +2 "$file1")
    mapfile -t lines2 < <(tail -n +2 "$file2")

    # Check if number of lines match
    if [[ ${#lines1[@]} -ne ${#lines2[@]} ]]; then
        return 1
    fi

    # Check each line
    for i in "${!lines1[@]}"; do
        local line1="${lines1[i]}"
        local line2="${lines2[i]}"

        # Split into words by tab
        IFS=$'\t' read -ra words1 <<< "$line1"
        IFS=$'\t' read -ra words2 <<< "$line2"

        # Check number of words
        if [[ ${#words1[@]} -ne 2 || ${#words2[@]} -ne 2 ]]; then
            return 1
        fi

        # Check first word
        if [[ "${words1[0]}" != "${words2[0]}" ]]; then
            return 1
        fi

        # Check second word: split by comma, sort, and compare
        IFS=',' read -ra nums1 <<< "${words1[1]}"
        IFS=',' read -ra nums2 <<< "${words2[1]}"

        # Sort and join
        local sorted1 sorted2
        IFS=$'\n' sorted1=($(sort <<<"${nums1[*]}"))
        IFS=$'\n' sorted2=($(sort <<<"${nums2[*]}"))

        if [[ "${sorted1[*]}" != "${sorted2[*]}" ]]; then
            return 1
        fi
    done

    return 0
}
