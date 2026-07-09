#!/usr/bin/env python3

import re
import io
import os
import csv
import json
import time
import argparse
import zipfile
import urllib.request
import urllib.error

import pandas as pd


OX_REGEX = re.compile(r"\bOX=(\d+)\b")
CORE_RANKS = {"domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Enrich one or more BLASTP TSV files with NCBI taxonomy"
    )
    parser.add_argument(
        "-i", "--input",
        nargs="+",
        help="One or more input blastp.tsv files"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV file (single input), output directory (multiple inputs), or lineage TSV in mode 1"
    )
    parser.add_argument(
        "--taxonomy",
        help="Intermediate taxonomy TSV file (required in mode 2)"
    )
    parser.add_argument(
        "--mode",
        choices=["all", "1", "2"],
        default="all",
        help="all: run everything; 1: build taxonomy only; 2: join using existing taxonomy"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100,
        help="NCBI taxonomy batch size"
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=1.0,
        help="Sleep time between NCBI batches"
    )
    return parser.parse_args()


def validate_args(args):
    if args.mode in {"all", "1", "2"} and not args.input:
        raise ValueError("--input is required")

    if args.mode == "2" and not args.taxonomy:
        raise ValueError("--taxonomy is required in mode 2")


def chunk_list(values, size):
    for i in range(0, len(values), size):
        yield values[i:i + size]


def parse_stitle(stitle):
    if pd.isna(stitle):
        return {
            "protein_name": None,
            "OS": None,
            "OX": None,
            "GN": None,
            "PE": None,
            "SV": None,
        }

    stitle = str(stitle).strip()

    pattern = re.compile(
        r"^(?P<protein_name>.*?)"
        r"(?:\s+OS=(?P<OS>.*?))?"
        r"(?:\s+OX=(?P<OX>\d+))?"
        r"(?:\s+GN=(?P<GN>.*?))?"
        r"(?:\s+PE=(?P<PE>\d+))?"
        r"(?:\s+SV=(?P<SV>\d+))?$"
    )

    m = pattern.match(stitle)
    if not m:
        return {
            "protein_name": stitle,
            "OS": None,
            "OX": None,
            "GN": None,
            "PE": None,
            "SV": None,
        }

    return m.groupdict()


def normalize_column(col):
    col = col.strip().lower()
    col = col.replace("/", "_")
    col = re.sub(r"\s+", "_", col)
    col = re.sub(r"[^a-z0-9_]", "", col)

    if col.startswith("domain_realm"):
        col = col.replace("domain_realm", "domain")

    return col


def is_core_rank_column(col):
    if col in {"OX_taxid", "ncbi_taxid", "ncbi_tax_name", "rank"}:
        return True

    if not (col.endswith("_name") or col.endswith("_taxid")):
        return False

    rank_part = col.rsplit("_", 1)[0]
    return any(rank_part == core or rank_part.endswith(core) for core in CORE_RANKS)


def standardize_ncbi_taxonomy_columns(df_tax):
    df_tax.columns = [normalize_column(c) for c in df_tax.columns]

    rename_map = {
        "query": "OX_taxid",
        "taxid": "ncbi_taxid",
        "tax_name": "ncbi_tax_name",
    }
    df_tax = df_tax.rename(columns={k: v for k, v in rename_map.items() if k in df_tax.columns})

    kept_cols = [c for c in df_tax.columns if is_core_rank_column(c)]
    df_tax = df_tax[kept_cols]

    if "OX_taxid" in df_tax.columns:
        df_tax["OX_taxid"] = df_tax["OX_taxid"].astype(str)

    return df_tax


def fetch_ncbi_taxonomy_batch(taxids, max_retries=4, sleep_time=3):
    url = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/download?filename=ncbi_dataset.zip"

    payload = {
        "tax_ids": [int(t) for t in taxids],
        "aux_reports": ["TAXONOMY_SUMMARY"]
    }
    payload_bytes = json.dumps(payload).encode("utf-8")

    headers = {
        "accept": "application/zip",
        "content-type": "application/json",
        "User-Agent": "blast-taxonomy-parser/1.0"
    }

    for attempt in range(max_retries):
        req = urllib.request.Request(url, data=payload_bytes, headers=headers, method="POST")

        try:
            with urllib.request.urlopen(req, timeout=120) as response:
                zip_bytes = response.read()

            with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
                target = "ncbi_dataset/data/taxonomy_summary.tsv"
                if target not in zf.namelist():
                    raise RuntimeError(f"{target} not found in NCBI zip archive")

                with zf.open(target) as fh:
                    df_tax = pd.read_csv(fh, sep="\t", dtype=str)

            return df_tax

        except urllib.error.HTTPError as e:
            print(f"Warning: HTTP {e.code} for batch of {len(taxids)} taxids")
        except Exception as e:
            print(f"Warning: failed NCBI batch of {len(taxids)} taxids: {e}")

        if attempt < max_retries - 1:
            wait = sleep_time * (attempt + 1)
            print(f"Retrying batch in {wait}s")
            time.sleep(wait)

    return None


def extract_taxids_from_file(filepath):
    taxids = set()

    with open(filepath, "r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)

        try:
            stitle_idx = header.index("stitle")
        except ValueError as e:
            raise RuntimeError(f"'stitle' column not found in {filepath}") from e

        for row in reader:
            if stitle_idx >= len(row):
                continue
            stitle = row[stitle_idx]
            m = OX_REGEX.search(stitle)
            if m:
                taxids.add(m.group(1))

    return taxids


def collect_all_taxids(filepaths):
    all_taxids = set()
    for fp in filepaths:
        print(f"Scanning taxids in {fp}")
        all_taxids.update(extract_taxids_from_file(fp))
    return sorted(all_taxids)


def build_taxonomy_df(all_taxids, batch_size, sleep_time):
    taxonomy_dfs = []
    failed_taxids = []

    for i, batch in enumerate(chunk_list(all_taxids, batch_size), start=1):
        print(f"Processing NCBI batch {i} ({len(batch)} taxids)")
        df_tax = fetch_ncbi_taxonomy_batch(batch)

        if df_tax is None or df_tax.empty:
            print(f"Warning: failed batch {i}")
            failed_taxids.extend(batch)
        else:
            taxonomy_dfs.append(df_tax)

        time.sleep(sleep_time)

    if taxonomy_dfs:
        taxonomy_df = pd.concat(taxonomy_dfs, ignore_index=True)
        taxonomy_df = standardize_ncbi_taxonomy_columns(taxonomy_df)
        taxonomy_df = taxonomy_df.drop_duplicates(subset=["OX_taxid"])
    else:
        taxonomy_df = pd.DataFrame(columns=["OX_taxid"])

    return taxonomy_df, sorted(set(failed_taxids))


def save_taxonomy_df(taxonomy_df, output_path):
    print(f"Writing taxonomy table: {output_path}")
    taxonomy_df.to_csv(output_path, sep="\t", index=False)


def load_taxonomy_df(taxonomy_path):
    print(f"Loading taxonomy table: {taxonomy_path}")
    df = pd.read_csv(taxonomy_path, sep="\t", dtype=str)
    if "OX_taxid" not in df.columns:
        raise RuntimeError(f"'OX_taxid' column not found in taxonomy file {taxonomy_path}")
    df["OX_taxid"] = df["OX_taxid"].astype(str)
    return df


def output_path_for_input(input_path, output_arg, multiple_inputs):
    if not multiple_inputs:
        return output_arg

    base = os.path.basename(input_path)
    if base.endswith(".tsv"):
        base = base[:-4] + "_with_tax.tsv"
    else:
        base = base + "_with_tax.tsv"

    return os.path.join(output_arg, base)


def process_one_file(input_path, output_path, taxonomy_df):
    print(f"Loading dataframe: {input_path}")
    df = pd.read_csv(input_path, sep="\t", dtype=str)

    parsed = df["stitle"].apply(parse_stitle).apply(pd.Series)
    df = pd.concat([df, parsed], axis=1)

    df["taxid"] = df["OX"].astype("string")

    if not taxonomy_df.empty:
        df = df.merge(
            taxonomy_df,
            how="left",
            left_on="taxid",
            right_on="OX_taxid"
        )

    print(f"Writing: {output_path}")
    df.to_csv(output_path, sep="\t", index=False)


def run_mode_1(args):
    all_taxids = collect_all_taxids(args.input)
    print(f"Unique taxids found across all files: {len(all_taxids)}")

    taxonomy_df, failed_taxids = build_taxonomy_df(
        all_taxids=all_taxids,
        batch_size=args.batch_size,
        sleep_time=args.sleep
    )

    save_taxonomy_df(taxonomy_df, args.output)

    if failed_taxids:
        failed_path = args.output + ".failed_taxids.tsv"
        pd.DataFrame({"taxid": failed_taxids}).to_csv(failed_path, sep="\t", index=False)
        print(f"Failed taxids saved to {failed_path}")


def run_mode_2(args):
    multiple_inputs = len(args.input) > 1
    taxonomy_df = load_taxonomy_df(args.taxonomy)

    if multiple_inputs:
        os.makedirs(args.output, exist_ok=True)

    for input_path in args.input:
        out_path = output_path_for_input(input_path, args.output, multiple_inputs)
        # check if output already exists
        if os.path.exists(out_path):
            print(f"Output already exists, skipping: {out_path}")
            continue
        process_one_file(input_path, out_path, taxonomy_df)


def run_mode_all(args):
    multiple_inputs = len(args.input) > 1

    if multiple_inputs:
        os.makedirs(args.output, exist_ok=True)

    all_taxids = collect_all_taxids(args.input)
    print(f"Unique taxids found across all files: {len(all_taxids)}")

    taxonomy_df, failed_taxids = build_taxonomy_df(
        all_taxids=all_taxids,
        batch_size=args.batch_size,
        sleep_time=args.sleep
    )

    lineage_path = os.path.join(args.output, "all_lineage.tsv") if multiple_inputs else args.output + ".all_lineage.tsv"
    save_taxonomy_df(taxonomy_df, lineage_path)

    if failed_taxids:
        failed_path = os.path.join(args.output, "failed_taxids.tsv") if multiple_inputs else args.output + ".failed_taxids.tsv"
        pd.DataFrame({"taxid": failed_taxids}).to_csv(failed_path, sep="\t", index=False)
        print(f"Failed taxids saved to {failed_path}")

    for input_path in args.input:
        out_path = output_path_for_input(input_path, args.output, multiple_inputs)
        process_one_file(input_path, out_path, taxonomy_df)


def main():
    args = parse_args()
    validate_args(args)

    if args.mode == "1":
        run_mode_1(args)
    elif args.mode == "2":
        run_mode_2(args)
    else:
        run_mode_all(args)


if __name__ == "__main__":
    main()