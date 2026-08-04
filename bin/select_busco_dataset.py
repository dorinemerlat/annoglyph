#!/usr/bin/env python3
"""
Select the most specific BUSCO dataset matching a taxonomy lineage.

The taxonomy file must be a tab-separated file with the columns:
rank, taxid, and name.

Usage:
    select_busco_dataset.py lineage.tsv
"""

import argparse
import re
import subprocess
from pathlib import Path


DATASET_PATTERN = re.compile(r"-\s+(\S+)\s+\[\d+\]")


def parse_args():
    parser = argparse.ArgumentParser(description="Select the most specific BUSCO dataset for a taxonomy lineage.")
    parser.add_argument("taxonomy", type=Path, help="Taxonomy lineage TSV file.")
    return parser.parse_args()


def list_busco_datasets():
    """Return the output of busco --list-datasets."""
    result = subprocess.run(
        ["busco", "--list-datasets"],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.splitlines()


def read_taxonomy_lineage(taxonomy_file):
    """Read lineage names from the most specific taxon to the root."""
    lineage = []

    with taxonomy_file.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\r\n")

            if not line:
                continue

            columns = line.split("\t")

            if line_number == 1 and columns[0] == "rank":
                continue

            if len(columns) < 3:
                raise ValueError(
                    f"Invalid taxonomy line {line_number}: expected at least 3 tab-separated columns."
                )

            lineage.append(columns[2].strip())

    return reversed(lineage)


def find_matching_dataset(lineage_names, dataset_lines):
    """Return the first BUSCO dataset matching the taxonomy lineage."""
    for taxon_name in lineage_names:
        for line in dataset_lines:
            if taxon_name.lower() not in line.lower():
                continue

            match = DATASET_PATTERN.search(line)

            if match:
                return match.group(1)

    return None


def main():
    args = parse_args()

    if not args.taxonomy.is_file():
        raise SystemExit(f"Taxonomy file does not exist: {args.taxonomy}")

    dataset_lines = list_busco_datasets()
    lineage_names = read_taxonomy_lineage(args.taxonomy)
    dataset = find_matching_dataset(lineage_names, dataset_lines)

    if not dataset:
        raise SystemExit("No matching BUSCO dataset found.")

    print(dataset)


if __name__ == "__main__":
    main()