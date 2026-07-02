#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $(basename "$0") output.tsv input1.tsv [input2.tsv ...]" >&2
    exit 1
fi

output="$1"
shift

awk '
FNR == 1 && NR != 1 { next }
{ print }
' "$@" > "$output"

if [[ ! -s "$output" ]]; then
    echo "ERROR: output file is empty: $output" >&2
    exit 1
fi
