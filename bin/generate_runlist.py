#!/usr/bin/env python3

import sys
import argparse
import requests
import re
import time
import xml.etree.ElementTree as ET


HEADER = "@Run_acc\ttotal_spots\ttotal_bases\tavg_len\tbool:paired\tcolor_space\n"


# ---------------------------------------------------------
# Small helper: HTTP GET with basic robustness
# ---------------------------------------------------------
def http_get(url: str, timeout: int = 60) -> str:
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return r.text


# ---------------------------------------------------------
# 1️⃣ VARUS-like retrieveID() — but DO NOT exit if not found
# ---------------------------------------------------------
def esearch_sra(genus: str, species: str):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=sra&term={genus}+{species}%5borgn%5d+AND+biomol_rna%5bProp%5d"
        "&usehistory=y"
    )

    try:
        xml = http_get(url)
    except Exception as e:
        print(f"ERROR: esearch request failed: {e}", file=sys.stderr)
        # No crash: return "no results"
        return None, 0

    root = ET.fromstring(xml)
    count = int(root.findtext(".//Count", default="0"))
    webenv = root.findtext(".//WebEnv")

    if count == 0 or not webenv:
        print("Phrase not found or no RNA datasets.", file=sys.stderr)
        return None, 0

    print(f"Server has {count} data sets", file=sys.stderr)
    return webenv, count


# ---------------------------------------------------------
# 2️⃣ VARUS-like getRuns() (URL style + regex parsing)
# ---------------------------------------------------------
def get_runs(webenv: str, retmax: int = 100, retstart: int = 0, only_paired: bool = False):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        f"?db=sra&WebEnv={webenv}&query_key=1"
        f"&retmax={retmax}&retstart={retstart}"
    )

    xml = http_get(url)
    lines = xml.splitlines()

    runs = []
    paired = 0
    colorspace = 0

    for line in lines:
        if re.search(r"LAYOUT.{1,10}PAIRED", line):
            paired = 1

        if re.search(r"Instrument ABI_SOLID", line):
            colorspace = 1

        m = re.search(r'Run acc="(.*?)" total_spots="(\d*?)" total_bases="(\d*?)"', line)
        if m:
            run_acc = m.group(1)
            total_spots = m.group(2)
            total_bases = m.group(3)

            if total_spots and total_bases:
                if (only_paired and paired) or (not only_paired):
                    avglen = int(100 * int(total_bases) / int(total_spots)) / 100
                    runs.append([run_acc, total_spots, total_bases, avglen, paired, colorspace])

            paired = 0
            colorspace = 0

    print(f"Found {len(runs)} runs in batch (retstart={retstart}, retmax={retmax})", file=sys.stderr)
    return runs


# ---------------------------------------------------------
# 3️⃣ Metadata for manual SRA IDs (kept close to VARUS parsing)
# ---------------------------------------------------------
def get_metadata_from_id(sra_id: str):
    # esearch to UID
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=sra&term={sra_id}&retmode=xml"
    )

    try:
        xml = http_get(url)
    except Exception:
        return None

    root = ET.fromstring(xml)
    uid = root.findtext(".//Id")
    if not uid:
        return None

    # esummary by UID (URL built manually)
    summary_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        f"?db=sra&id={uid}"
    )

    try:
        text = http_get(summary_url)
    except Exception:
        return None

    lines = text.splitlines()
    paired = 0
    colorspace = 0

    for line in lines:
        if re.search(r"LAYOUT.{1,10}PAIRED", line):
            paired = 1
        if re.search(r"Instrument ABI_SOLID", line):
            colorspace = 1

        m = re.search(r'Run acc="(.*?)" total_spots="(\d*?)" total_bases="(\d*?)"', line)
        if m:
            run_acc = m.group(1)
            total_spots = m.group(2)
            total_bases = m.group(3)

            if total_spots and total_bases:
                avglen = int(100 * int(total_bases) / int(total_spots)) / 100
                return [run_acc, total_spots, total_bases, avglen, paired, colorspace]

    return None


# ---------------------------------------------------------
# 4️⃣ Main
# ---------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Generate VARUS-style Runlist")

    parser.add_argument("-g", "--genus", required=True)
    parser.add_argument("-s", "--specie", required=True)
    parser.add_argument("-i", "--sra-ids", help="Comma-separated SRA IDs")
    parser.add_argument("--only-paired", action="store_true")
    parser.add_argument("--retmax", type=int, default=100, help="Batch size for esummary (default: 100)")
    parser.add_argument("--retstart", type=int, default=0, help="Start offset for esummary (default: 0)")
    parser.add_argument("--all", action="store_true", help="Retrieve all available runs (paginate)")
    parser.add_argument("--sleep", type=float, default=0.3, help="Sleep between requests (default: 0.3s)")
    parser.add_argument("-o", "--output", help="Output file (default stdout)")

    args = parser.parse_args()

    out = open(args.output, "w") if args.output else sys.stdout
    out.write(HEADER)

    # ---- PART 1: species search (VARUS logic) but never crash pipeline
    webenv, count = esearch_sra(args.genus, args.specie)

    if webenv:
        if args.all:
            # paginate through all results
            retstart = args.retstart
            while retstart < count:
                batch = get_runs(
                    webenv,
                    retmax=args.retmax,
                    retstart=retstart,
                    only_paired=args.only_paired,
                )
                for r in batch:
                    out.write("\t".join(map(str, r)) + "\n")
                retstart += args.retmax
                time.sleep(args.sleep)
        else:
            runs = get_runs(
                webenv,
                retmax=args.retmax,
                retstart=args.retstart,
                only_paired=args.only_paired,
            )
            for r in runs:
                out.write("\t".join(map(str, r)) + "\n")

    # ---- PART 2: manual SRA IDs (always add if provided)
    if args.sra_ids:
        extra_ids = [x.strip() for x in args.sra_ids.split(",") if x.strip()]
        for sra_id in extra_ids:
            meta = get_metadata_from_id(sra_id)
            if meta:
                out.write("\t".join(map(str, meta)) + "\n")
            else:
                out.write(f"{sra_id}\t0\t0\t0\t0\t0\n")
            time.sleep(args.sleep)

    if args.output:
        out.close()


if __name__ == "__main__":
    main()