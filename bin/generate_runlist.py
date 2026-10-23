#!/usr/bin/env python3
"""Generate SRA metadata and VARUS/STAR runlists."""

import argparse
import csv
import sys
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import requests


EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

RUNLIST_COLUMNS = [
    "@Run_acc",
    "total_spots",
    "total_bases",
    "avg_len",
    "bool:paired",
    "color_space",
]

METADATA_COLUMNS = [
    "Run_acc",
    "total_spots",
    "total_bases",
    "avg_len",
    "bool:paired",
    "color_space",
    "experiment_accession",
    "study_accession",
    "study_title",
    "sample_accession",
    "sample_name",
    "experiment_title",
    "library_name",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "instrument_model",
    "platform",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Generate SRA metadata and VARUS/STAR runlists.")
    parser.add_argument("-g", "--genus", required=True, help="Genus name.")
    parser.add_argument("-s", "--species", required=True, help="Species epithet.")
    parser.add_argument("-i", "--sra-ids", help="Comma-separated SRA run accessions.")
    parser.add_argument("--only-paired", action="store_true", help="Keep only paired-end runs in the filtered STAR runlist.")
    parser.add_argument("--limit", type=int, default=0, help="Maximum number of runs in the filtered STAR runlist; 0 keeps all eligible runs.")
    parser.add_argument("--max-bases", type=int, default=0, help="Maximum number of bases per run in the filtered STAR runlist; 0 disables the filter.")
    parser.add_argument("--retmax", type=int, default=100, help="EUtils batch size.")
    parser.add_argument("--sleep", type=float, default=0.34, help="Delay between EUtils requests.")
    parser.add_argument("--runlist-output", type=Path, default=Path("Runlist.txt"), help="Complete VARUS-compatible runlist.")
    parser.add_argument("--filtered-runlist-output", type=Path, default=Path("Runlist_filtered.txt"), help="Filtered runlist used for STAR.")
    parser.add_argument("--metadata-output", type=Path, default=Path("run_metadata.tsv"), help="Complete metadata output.")
    return parser.parse_args()


def http_get(endpoint, params, timeout=60):
    """Retrieve an NCBI EUtils response."""
    response = requests.get(f"{EUTILS_URL}/{endpoint}", params=params, timeout=timeout)
    response.raise_for_status()
    return response.text


def text_or_empty(element, path):
    """Return stripped XML text or an empty string."""
    if element is None:
        return ""

    value = element.findtext(path, default="")
    return value.strip() if value else ""


def first_text(element, paths):
    """Return the first non-empty XML value."""
    for path in paths:
        value = text_or_empty(element, path)
        if value:
            return value

    return ""


def esearch_sra(genus, species):
    """Search RNA-seq experiments for one species."""
    term = f"{genus} {species}[Organism] AND biomol rna[Properties]"
    xml = http_get(
        "esearch.fcgi",
        {
            "db": "sra",
            "term": term,
            "usehistory": "y",
            "retmode": "xml",
        },
    )

    root = ET.fromstring(xml)
    count = int(root.findtext(".//Count", default="0"))
    webenv = root.findtext(".//WebEnv", default="")
    query_key = root.findtext(".//QueryKey", default="")

    print(f"Found {count} SRA experiments for {genus} {species}.", file=sys.stderr)
    return webenv, query_key, count


def esearch_accession(accession):
    """Resolve one SRA accession to an Entrez UID."""
    xml = http_get(
        "esearch.fcgi",
        {
            "db": "sra",
            "term": f"{accession}[Accession]",
            "retmode": "xml",
        },
    )

    root = ET.fromstring(xml)
    return root.findtext(".//Id", default="")


def parse_fragment(text):
    """Parse an XML fragment stored inside an ESummary field."""
    if not text:
        return None

    try:
        return ET.fromstring(text)
    except ET.ParseError:
        try:
            return ET.fromstring(f"<root>{text}</root>")
        except ET.ParseError:
            return None


def get_docsum_item(docsum, name):
    """Return one named ESummary item."""
    item = docsum.find(f"./Item[@Name='{name}']")
    return item.text if item is not None and item.text else ""


def parse_docsum(docsum):
    """Extract runs and experiment metadata from one SRA DocSum."""
    exp_xml = parse_fragment(get_docsum_item(docsum, "ExpXml"))
    runs_xml = parse_fragment(get_docsum_item(docsum, "Runs"))

    metadata = {
        "experiment_accession": "",
        "study_accession": "",
        "study_title": "",
        "sample_accession": "",
        "sample_name": "",
        "experiment_title": "",
        "library_name": "",
        "library_strategy": "",
        "library_source": "",
        "library_selection": "",
        "library_layout": "",
        "instrument_model": "",
        "platform": "",
    }

    if exp_xml is not None:
        experiment = exp_xml.find(".//Experiment")
        if experiment is None:
            experiment = exp_xml.find(".//EXPERIMENT")
        study = exp_xml.find(".//Study")
        if study is None:
            study = exp_xml.find(".//STUDY")
        sample = exp_xml.find(".//Sample")
        if sample is None:
            sample = exp_xml.find(".//SAMPLE")

        if experiment is not None:
            metadata["experiment_accession"] = experiment.attrib.get("acc", "")

        if study is not None:
            metadata["study_accession"] = study.attrib.get("acc", "")
            metadata["study_title"] = study.attrib.get("name", "")

        if sample is not None:
            metadata["sample_accession"] = sample.attrib.get("acc", "")
            metadata["sample_name"] = sample.attrib.get("name", "")

        metadata["experiment_title"] = first_text(
            exp_xml,
            [
                ".//TITLE",
                ".//Title",
                ".//EXPERIMENT_TITLE",
            ],
        )
        metadata["library_name"] = first_text(exp_xml, [".//LIBRARY_NAME", ".//Library_name"])
        metadata["library_strategy"] = first_text(exp_xml, [".//LIBRARY_STRATEGY", ".//Library_strategy"])
        metadata["library_source"] = first_text(exp_xml, [".//LIBRARY_SOURCE", ".//Library_source"])
        metadata["library_selection"] = first_text(exp_xml, [".//LIBRARY_SELECTION", ".//Library_selection"])
        metadata["instrument_model"] = first_text(exp_xml, [".//INSTRUMENT_MODEL", ".//Instrument_model"])
        summary_platform = exp_xml.find(".//Platform")
        if summary_platform is not None:
            metadata["instrument_model"] = metadata["instrument_model"] or summary_platform.get("instrument_model", "")
            metadata["platform"] = (summary_platform.text or "").strip()

        layout = exp_xml.find(".//LIBRARY_LAYOUT")
        if layout is not None and list(layout):
            metadata["library_layout"] = list(layout)[0].tag.upper()

        platform = exp_xml.find(".//PLATFORM")
        if platform is not None and list(platform):
            metadata["platform"] = list(platform)[0].tag.upper()

    rows = []

    if runs_xml is None:
        return rows

    # ESummary stores a single run as the root element, not inside a wrapper.
    for run in (element for element in runs_xml.iter() if element.tag in ("Run", "RUN")):
        accession = run.attrib.get("acc", "")
        total_spots = int(run.attrib.get("total_spots") or 0)
        total_bases = int(run.attrib.get("total_bases") or 0)

        if not accession:
            continue

        avg_len = round(total_bases / total_spots, 2) if total_spots else 0
        paired = int(metadata["library_layout"].upper() == "PAIRED")
        color_space = int("SOLID" in metadata["instrument_model"].upper())

        rows.append(
            {
                "Run_acc": accession,
                "total_spots": total_spots,
                "total_bases": total_bases,
                "avg_len": avg_len,
                "bool:paired": paired,
                "color_space": color_space,
                **metadata,
            }
        )

    return rows


def fetch_docsum_batch(params):
    """Retrieve and parse one ESummary batch."""
    xml = http_get("esummary.fcgi", params)
    root = ET.fromstring(xml)
    rows = []

    for docsum in root.findall(".//DocSum"):
        rows.extend(parse_docsum(docsum))

    return rows


def fetch_species_runs(webenv, query_key, count, retmax, sleep):
    """Retrieve all runs matching the species search."""
    rows = []

    for retstart in range(0, count, retmax):
        batch = fetch_docsum_batch(
            {
                "db": "sra",
                "WebEnv": webenv,
                "query_key": query_key,
                "retstart": retstart,
                "retmax": retmax,
                "retmode": "xml",
            }
        )
        rows.extend(batch)
        print(f"Retrieved {len(batch)} runs from batch starting at {retstart}.", file=sys.stderr)
        time.sleep(sleep)

    return rows


def fetch_accession_runs(accession):
    """Retrieve metadata for one manually supplied accession."""
    uid = esearch_accession(accession)

    if not uid:
        print(f"WARNING: SRA accession not found: {accession}", file=sys.stderr)
        return []

    return fetch_docsum_batch(
        {
            "db": "sra",
            "id": uid,
            "retmode": "xml",
        }
    )


def deduplicate_runs(rows):
    """Keep one metadata row per run accession."""
    by_accession = {}

    for row in rows:
        accession = row["Run_acc"]

        if accession not in by_accession:
            by_accession[accession] = row
            continue

        current = by_accession[accession]
        current_score = sum(bool(current.get(column)) for column in METADATA_COLUMNS)
        new_score = sum(bool(row.get(column)) for column in METADATA_COLUMNS)

        if new_score > current_score:
            by_accession[accession] = row

    return list(by_accession.values())


def select_runs(rows, limit, max_bases, only_paired):
    """Select a diverse and technically suitable subset for STAR."""
    eligible = [
        row
        for row in rows
        if not row["color_space"]
        and (max_bases <= 0 or row["total_bases"] <= max_bases)
        and (not only_paired or row["bool:paired"])
    ]

    eligible.sort(
        key=lambda row: (
            -row["bool:paired"],
            row["study_accession"],
            row["sample_accession"],
            -row["total_bases"],
            row["Run_acc"],
        )
    )

    if limit <= 0 or len(eligible) <= limit:
        return eligible

    grouped = defaultdict(list)

    for row in eligible:
        group = (
            row["study_accession"] or "unknown_study",
            row["sample_accession"] or row["Run_acc"],
        )
        grouped[group].append(row)

    selected = []
    groups = list(grouped.values())

    while groups and len(selected) < limit:
        remaining_groups = []

        for group in groups:
            if group and len(selected) < limit:
                selected.append(group.pop(0))

            if group:
                remaining_groups.append(group)

        groups = remaining_groups

    return selected


def write_tsv(path, rows, columns):
    """Write rows as a TSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_runlist(path, rows):
    """Write a VARUS-compatible runlist."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RUNLIST_COLUMNS, delimiter="\t")
        writer.writeheader()

        for row in rows:
            writer.writerow(
                {
                    "@Run_acc": row["Run_acc"],
                    "total_spots": row["total_spots"],
                    "total_bases": row["total_bases"],
                    "avg_len": row["avg_len"],
                    "bool:paired": row["bool:paired"],
                    "color_space": row["color_space"],
                }
            )


def main():
    args = parse_args()

    webenv, query_key, count = esearch_sra(args.genus, args.species)
    rows = []

    if webenv and query_key and count:
        rows.extend(
            fetch_species_runs(
                webenv,
                query_key,
                count,
                args.retmax,
                args.sleep,
            )
        )

    if args.sra_ids:
        accessions = [value.strip() for value in args.sra_ids.split(",") if value.strip()]

        for accession in accessions:
            rows.extend(fetch_accession_runs(accession))
            time.sleep(args.sleep)

    rows = deduplicate_runs(rows)
    rows.sort(key=lambda row: row["Run_acc"])

    selected_rows = select_runs(
        rows,
        limit=args.limit,
        max_bases=args.max_bases,
        only_paired=args.only_paired,
    )

    oversized_runs = [
        row
        for row in rows
        if args.max_bases > 0 and row["total_bases"] > args.max_bases
    ]

    if oversized_runs:
        print(
            f"Excluded {len(oversized_runs)} runs exceeding {args.max_bases:,} bases.",
            file=sys.stderr,
        )

    write_tsv(args.metadata_output, rows, METADATA_COLUMNS)
    write_runlist(args.runlist_output, rows)
    write_runlist(args.filtered_runlist_output, selected_rows)

    print(f"Wrote metadata for {len(rows)} runs to {args.metadata_output}.", file=sys.stderr)
    print(f"Wrote {len(rows)} runs to the complete VARUS runlist: {args.runlist_output}.", file=sys.stderr)
    print(f"Selected {len(selected_rows)} runs for STAR: {args.filtered_runlist_output}.", file=sys.stderr)


if __name__ == "__main__":
    main()