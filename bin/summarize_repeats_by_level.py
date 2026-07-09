#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys
import gffpandas.gffpandas as gffpd
import pandas as pd
import os
import subprocess

def parse_args():
    p = argparse.ArgumentParser(
        description="Convert RepeatMasker .out to GFF with 9 columns"
    )
    p.add_argument('-i', '--input', required=True, help='Gff file')
    p.add_argument('-s', '--size', required=True, type=float, help='Genome size')
    return p.parse_args()


def parse_attributes(s):
    if pd.isna(s) or s == "":
        return {}
    attributes = {}
    for item in str(s).split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            attributes[k] = v
    return attributes


def add_genome_id(path, df):
    genome_id = os.path.basename(path).split(".")[0].split("_")[0]
    df.insert(0, "genome_id", genome_id)
    return df, genome_id

def parse_gff(path):
    """Parse a GFF file into a DataFrame and expand attributes into columns."""
    gff_obj = gffpd.read_gff3(path)

    gff = gff_obj.df.copy()

    # Parse attributes into dict
    attr_dicts = gff["attributes"].apply(parse_attributes)

    # Expand dicts into columns
    attr_df = pd.json_normalize(attr_dicts)

    # Merge expanded attributes with the main dataframe
    gff = gff.reset_index(drop=True)

    attr_dicts = gff["attributes"].apply(parse_attributes)
    attr_df = pd.json_normalize(attr_dicts).reset_index(drop=True)

    gff = pd.concat([gff.drop(columns=["attributes"]), attr_df], axis=1)

    return gff


def calculate_coverage(gff):
    bed_cols = ["seq_id", "start", "end"]
    bed = gff[bed_cols]

    # write bed to temp file
    with open("temp.bed", "w") as bed_file:
        bed.to_csv(bed_file, sep="\t", index=False, header=False)
    # sort and merge bed using bedtools
    merged_bed_path = "merged_temp.bed"
    cmd = "awk '{ if ($2 > $3) { t=$2; $2=$3; $3=t } print }' OFS='\t' temp.bed| sort -k1,1 -k2,2n| /biolo/bedtools/2.31.1/bin/bedtools merge -i -"
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
    covered_bp = 0

    for line in p.stdout.splitlines():
        chrom, start, end = line.split("\t")
        covered_bp += int(end) - int(start) + 1

    return covered_bp


def main():
    print("Starting")
    args = parse_args()
    in_path = Path(args.input)
    if not in_path.exists():
        print(f"Error: input .out not found: {in_path}", file=sys.stderr)
        sys.exit(2)

    gff = parse_gff(in_path)
    gff = gff.rename(columns={"type": "feature"})
    gff = gff.rename(columns=lambda c: c.replace("repeat_", ""))

    # filter to keep only relevant columns
    bed_cols = [
        "seq_id",
        "start",
        "end"
    ]

    levels = [
        "feature",
        "type",
        "class",
        "order",
        "superfamily",
        "family",
        "subfamily",
        "subfamily2"
            ]

    # --- Extra summary: all_only_repetitive_elements (exclude pseudogene + low_complexity) ---
    excluded_types = {"pseudogene", "low_complexity"}

    if "type" not in gff.columns:
        print("Warning: column 'type' not found; cannot compute all_only_repetitive_elements", file=sys.stderr)
    else:
        gf_rep_only = gff.loc[~gff["type"].isin(excluded_types), bed_cols].copy()

        print(f"Summarizing by: {bed_cols} (all_only_repetitive_elements)")

        results = []
        covered_bp = calculate_coverage(gf_rep_only)
        load = len(gf_rep_only)
        results.append({"covered_bp": covered_bp, "load": load})

        summary_df = pd.DataFrame(results)
        summary_df, genome_id = add_genome_id(in_path, summary_df)

        genome_size = float(args.size)
        summary_df["coverage"] = summary_df["covered_bp"] * 100.0 / genome_size

        out_file = f"{genome_id}_repeats_summary_by_all_only_repetitive_elements.tsv"
        summary_df.to_csv(out_file, sep="\t", index=False)
        print(f"Saved: {out_file}")


    for level in range(len(levels) + 1 ):
        cols_to_group = levels[: level ]
        kept_cols = bed_cols + cols_to_group
        print(f"Summarizing by: {kept_cols}")

        # filter columns
        gf_filtered = gff[bed_cols + cols_to_group]

        # extract element types
        sub_uniq = gff[cols_to_group].drop_duplicates()

        # Build groups (one group if no level column)
        results = []

        if len(cols_to_group) == 0:
            # single group = all rows
            covered_bp = calculate_coverage(gf_filtered)
            load = len(gf_filtered)
            results.append({"covered_bp": covered_bp, "load": load})

        else:
            for keys, subdf in gf_filtered.groupby(cols_to_group, dropna=False, sort=False):
                if not isinstance(keys, tuple):
                    keys = (keys,)

                row = {col: val for col, val in zip(cols_to_group, keys)}

                covered_bp = calculate_coverage(subdf)
                row["covered_bp"] = covered_bp

                # load = number of occurrences (rows) in this group
                row["load"] = len(subdf)

                results.append(row)

        summary_df = pd.DataFrame(results)
        summary_df, genome_id = add_genome_id(in_path, summary_df)

        # add coverage %
        genome_size = float(args.size)
        summary_df["coverage"] = summary_df["covered_bp"] * 100.0 / genome_size

        # nice deterministic output order
        if len(cols_to_group) > 0 and not summary_df.empty:
            summary_df = summary_df.sort_values(cols_to_group, ascending=False).reset_index(drop=True)

        # save dataframe
        level_name = "all" if len(cols_to_group) == 0 else "__".join(cols_to_group)
        out_file = f"{genome_id}_repeats_summary_by_{level_name}.tsv"

        summary_df.to_csv(out_file, sep="\t", index=False)
        print(f"Saved: {out_file}")

if __name__ == "__main__":
    main()
