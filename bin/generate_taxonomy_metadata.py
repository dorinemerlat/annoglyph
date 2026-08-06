#!/usr/bin/env python3

import argparse
import json
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

import pandas as pd


NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UNIPROT_URL = "https://rest.uniprot.org/taxonomy/{taxid}"

EXCLUDED_RANKS = {"no rank", "clade"}


def parse_args():
    parser = argparse.ArgumentParser(description="Generate taxonomy metadata.")
    parser.add_argument("-i", "--input", required=True, type=Path, help="Input TSV.")
    parser.add_argument("-o", "--output", default=Path("output"), type=Path, help="Output directory.")
    return parser.parse_args()


def fetch(url, accept):
    """Download a remote resource."""

    request = urllib.request.Request(
        url,
        headers={
            "Accept": accept,
            "User-Agent": "EXOGAP",
        },
    )

    with urllib.request.urlopen(request) as response:
        return response.read()


def fetch_taxonomy(taxid):
    """Return NCBI taxonomy information."""

    query = urllib.parse.urlencode(
        {
            "db": "taxonomy",
            "id": taxid,
            "retmode": "xml",
        }
    )

    root = ET.fromstring(
        fetch(f"{NCBI_EFETCH_URL}?{query}", "application/xml")
    )

    taxon = root.find("Taxon")

    lineage = []

    for node in taxon.findall("./LineageEx/Taxon") + [taxon]:

        rank = node.findtext("Rank")

        if rank in EXCLUDED_RANKS:
            continue

        lineage.append(
            {
                "rank": rank,
                "taxid": node.findtext("TaxId"),
                "name": node.findtext("ScientificName"),
            }
        )

    return {
        "taxid": taxon.findtext("TaxId"),
        "scientific_name": taxon.findtext("ScientificName"),
        "lineage": lineage,
    }


def fetch_mnemonic(taxid):
    """Return the UniProt mnemonic."""

    data = json.loads(
        fetch(
            UNIPROT_URL.format(taxid=taxid),
            "application/json",
        )
    )

    return str(data.get("mnemonic") or "NA")


def fallback_mnemonic(genome_id):
    """Generate a mnemonic from the genome identifier."""

    parts = genome_id.split("-")

    if len(parts) >= 2:
        return (parts[0][:3] + parts[1][:2]).upper()

    return re.sub(r"[^A-Za-z0-9]", "", genome_id)[:5].upper()


def assign_mnemonics(records):
    """Assign official and unique EXOGAP mnemonics."""

    bases = []

    for record in records:

        user = record["user_mnemonic"]
        official = record["official_mnemonic"]

        if user and user != "NA":
            base = user
            official = user

        elif official != "NA" and not official.startswith("9"):
            base = official

        else:
            base = fallback_mnemonic(record["id"])

        record["official_mnemonic"] = official
        record["_base"] = base
        bases.append(base)

    counts = Counter(bases)
    index = Counter()

    for record in records:

        base = record.pop("_base")

        if counts[base] == 1:
            record["exogap_mnemonic"] = base
        else:
            index[base] += 1
            record["exogap_mnemonic"] = f"{base}{index[base]}"

def build_tree(records):
    """Build a taxonomy tree."""

    tree = {}

    for record in records:

        node = tree

        for taxon in record["taxonomy"]["lineage"]:

            node = node.setdefault(
                taxon["taxid"],
                {
                    "name": taxon["name"],
                    "children": {},
                },
            )["children"]

    return tree


def write_newick(node):
    """Convert a taxonomy tree to Newick."""

    children = []

    for child in node.values():

        subtree = write_newick(child["children"])

        if subtree:
            children.append(f"({subtree}){child['name']}")
        else:
            children.append(child["name"])

    return ",".join(children)


def main():

    args = parse_args()
    args.output.mkdir(exist_ok=True)

    df = pd.read_csv(args.input, sep="\t", dtype=str).fillna("")

    taxonomy_cache = {}
    mnemonic_cache = {}

    records = []
    used_ranks = []

    for _, row in df.iterrows():

        taxid = row["taxid"]

        if taxid not in taxonomy_cache:
            taxonomy_cache[taxid] = fetch_taxonomy(taxid)

        if taxid not in mnemonic_cache:
            mnemonic_cache[taxid] = fetch_mnemonic(taxid)

        taxonomy = taxonomy_cache[taxid]

        for taxon in taxonomy["lineage"]:

            if taxon["rank"] not in used_ranks:
                used_ranks.append(taxon["rank"])

        records.append(
            {
                "id": row["id"],
                "name": row["name"],
                "taxonomy": taxonomy,
                "user_mnemonic": row.get("mnemonic", ""),
                "official_mnemonic": mnemonic_cache[taxid],
            }
        )

    assign_mnemonics(records)

    tree = build_tree(records)

    with open(args.output / "tree.nwk", "w") as out:
        out.write(f"({write_newick(tree)});\n")