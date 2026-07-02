#!/usr/bin/env bash
set -euo pipefail

taxonomy_file="$1"
lineage_file="lineage.txt"

busco --list-datasets > "$lineage_file"

while IFS=$'\t' read -r rank taxid name_rest; do
    [[ -z "${rank:-}" ]] && continue
    [[ "$rank" == "rank" ]] && continue

    group_name=$(echo "$name_rest" | tr -d '\r')
    line=$(grep -i "$group_name" "$lineage_file" | head -n 1 || true)
    match=$(echo "$line" | sed -E 's/.*- ([^ ]+) \[[0-9]+\].*/\1/')

    if [[ -n "$match" && "$match" != "$line" ]]; then
        echo "$match"
        exit 0
    fi
done < <(tac "$taxonomy_file")

echo "No matching BUSCO dataset found" >&2
exit 1