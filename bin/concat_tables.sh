#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: concat_table.sh output.tsv input1.tsv [input2.tsv ...]"
    exit 1
fi

output="$1"
shift

first_file="$1"

# Write header
head -n 1 "$first_file" > "$output"

# Append data from all files, skipping headers
for file in "$@"; do
    tail -n +2 "$file" >> "$output"
done

# Check output is not empty
[[ -s "$output" ]] || {
    echo "ERROR: Output file '$output' is empty."
    exit 1
}