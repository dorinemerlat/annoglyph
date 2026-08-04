#!/usr/bin/env python3
"""
Build a Newick taxonomy tree from species names and NCBI Taxonomy identifiers.

Usage:
    download_newick.py [--outfile FILE] species/taxid [...]

Examples:
    download_newick.py Homo_sapiens/9606 Mus_musculus/10090
    download_newick.py --outfile tree.nwk Homo_sapiens/9606 Mus_musculus/10090
"""

import argparse
import re
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=SyntaxWarning)

from ete3 import NCBITaxa


def parse_args():
    parser = argparse.ArgumentParser(description="Build a Newick tree from NCBI Taxonomy identifiers.")
    parser.add_argument("taxa", nargs="+", help="Taxa formatted as species/taxid.")
    parser.add_argument("--outfile", type=Path, help="Output Newick file. Defaults to stdout.")
    parser.add_argument("--database", type=Path, default=Path("taxa.sqlite"), help="ETE3 taxonomy database path.")
    parser.add_argument("--with-distances", action="store_true", help="Keep branch lengths in the Newick tree.")
    return parser.parse_args()


def parse_taxa(items):
    """Parse species/taxid entries into a taxid-to-name mapping."""
    taxa = {}

    for item in items:
        if "/" not in item:
            raise ValueError(f"Invalid taxon entry: {item}. Expected species/taxid.")

        species, taxid = item.rsplit("/", 1)

        if not taxid.isdigit():
            raise ValueError(f"Invalid TaxID in entry: {item}.")

        taxa.setdefault(int(taxid), []).append(species)

    return taxa


def rename_leaves(newick, name_map):
    """Replace leaf TaxIDs with species_name|taxid labels."""
    for taxid, names in name_map.items():
        label = f"{names[0]}|{taxid}"
        newick = re.sub(rf"(?<=[(,]){taxid}(?=[:,)])", label, newick)

    return newick


def rename_internal_nodes(newick, ncbi):
    """Replace internal-node TaxIDs with scientific_name|taxid labels."""
    taxids = {int(taxid) for taxid in re.findall(r"\)(\d+)(?=:|,|\)|;)", newick)}

    if not taxids:
        return newick

    names = ncbi.get_taxid_translator(taxids)

    for taxid, name in names.items():
        label = f"{name.replace(' ', '-')}|{taxid}"
        newick = re.sub(rf"(?<=\)){taxid}(?=:|,|\)|;)", label, newick)

    return newick


def build_newick(taxids, database):
    """Build a taxonomy tree using the local ETE3 taxonomy database."""
    ncbi = NCBITaxa(dbfile=str(database))
    tree = ncbi.get_topology(taxids)
    return tree.write(format=1, format_root_node=True), ncbi


def main():
    args = parse_args()
    name_map = parse_taxa(args.taxa)

    newick, ncbi = build_newick(list(name_map), args.database)
    newick = rename_leaves(newick, name_map)
    newick = rename_internal_nodes(newick, ncbi)

    if not args.with_distances:
        newick = re.sub(r":1(?:\.0+)?", "", newick)

    if args.outfile:
        args.outfile.write_text(f"{newick}\n", encoding="utf-8")
        print(f"Newick tree saved to {args.outfile}")
    else:
        print(newick)


if __name__ == "__main__":
    main()