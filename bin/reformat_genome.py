#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description="Rename FASTA sequence identifiers and write an ID mapping table.")
    parser.add_argument("--fasta", required=True, type=Path, help="Input FASTA file.")
    parser.add_argument("--prefix", required=True, help="Prefix used for the new sequence identifiers.")
    parser.add_argument("--output", required=True, type=Path, help="Output FASTA file.")
    parser.add_argument("--mapping", required=True, type=Path, help="Output TSV mapping file.")
    return parser.parse_args()


def main():
    args = parse_args()

    record_count = 0

    with args.output.open("w") as fasta_handle, args.mapping.open("w", newline="") as mapping_handle:
        writer = csv.writer(mapping_handle, delimiter="\t")
        writer.writerow(["previous_IDs", "EXOGAP_IDs"])

        for record_count, record in enumerate(SeqIO.parse(args.fasta, "fasta"), start=1):
            previous_id = record.id
            new_id = f"{args.prefix}_{record_count}"

            writer.writerow([previous_id, new_id])

            record.id = new_id
            record.name = new_id
            record.description = ""
            SeqIO.write(record, fasta_handle, "fasta")

    if record_count == 0:
        raise ValueError(f"No FASTA records found in {args.fasta}")


if __name__ == "__main__":
    main()