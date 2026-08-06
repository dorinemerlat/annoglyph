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

    root = ET.fromstring(fetch(f"{NCBI_EFETCH_URL}?{query}", "application/xml"))
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
    data = json.loads(fetch(UNIPROT_URL.format(taxid=taxid), "application/json"))
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
    indices = Counter()

    for record in records:
        base = record.pop("_base")

        if counts[base] == 1:
            record["exogap_mnemonic"] = base
        else:
            indices[base] += 1
            record["exogap_mnemonic"] = f"{base}{indices[base]}"


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


def write_newick(nodes):
    """Convert a taxonomy tree to Newick."""
    children = []

    for node in nodes.values():
        subtree = write_newick(node["children"])

        if subtree:
            children.append(f"({subtree}){node['name']}")
        else:
            children.append(node["name"])

    return ",".join(children)


def write_minimal_newick(nodes):
    """Convert a taxonomy tree to Newick and collapse single-child nodes."""
    children = []

    for node in nodes.values():
        node_children = node["children"]

        if not node_children:
            children.append(node["name"])
            continue

        subtree = write_minimal_newick(node_children)

        if len(node_children) == 1:
            children.append(subtree)
        else:
            children.append(f"({subtree}){node['name']}")

    return ",".join(children)


def main():
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    dataframe = pd.read_csv(args.input, sep="\t", dtype=str).fillna("")

    taxonomy_cache = {}
    mnemonic_cache = {}
    records = []
    used_ranks = []

    for _, row in dataframe.iterrows():
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

    rows = []
    taxonomy_json = {}

    for record in records:
        taxonomy = record["taxonomy"]

        row = {
            "id": record["id"],
            "name": record["name"],
            "taxonomy": taxonomy["scientific_name"],
            "taxid": taxonomy["taxid"],
            "official_mnemonic": record["official_mnemonic"],
            "exogap_mnemonic": record["exogap_mnemonic"],
        }

        lineage = {}

        for taxon in taxonomy["lineage"]:
            row[taxon["rank"]] = taxon["name"]
            lineage[taxon["rank"]] = {
                "name": taxon["name"],
                "taxid": taxon["taxid"],
            }

        rows.append(row)

        taxonomy_json[record["id"]] = {
            "id": record["id"],
            "name": record["name"],
            "taxonomy": taxonomy["scientific_name"],
            "taxid": taxonomy["taxid"],
            "mnemonic": {
                "official": record["official_mnemonic"],
                "exogap": record["exogap_mnemonic"],
            },
            "lineage": lineage,
        }

    columns = [
        "id",
        "name",
        "taxonomy",
        "taxid",
        "official_mnemonic",
        "exogap_mnemonic",
        *used_ranks,
    ]

    pd.DataFrame(rows).reindex(columns=columns).to_csv(
        args.output / "taxonomy.tsv",
        sep="\t",
        index=False,
    )

    with open(args.output / "taxonomy.json", "w", encoding="utf-8") as output:
        json.dump(taxonomy_json, output, indent=2, ensure_ascii=False)
        output.write("\n")

    tree = build_tree(records)

    with open(args.output / "tree.nwk", "w", encoding="utf-8") as output:
        output.write(f"({write_newick(tree)});\n")

    with open(args.output / "tree_minimal.nwk", "w", encoding="utf-8") as output:
        output.write(f"({write_minimal_newick(tree)});\n")


if __name__ == "__main__":
    main()