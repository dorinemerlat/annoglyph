#!/usr/bin/env python3
"""
Simple utility to download a Newick tree from a list of species/taxid,
with all nodes (feuilles et internes) étiquetés Name|taxID.

Usage:
    download_newick.py [--outfile FILE] species/taxid [...]

Examples:
    download_newick.py Homo_sapiens/9606 Mus_musculus/10090 Drosophila_melanogaster/7227
    download_newick.py --outfile tree.nwk Homo_sapiens/9606 Mus_musculus/10090
"""

import argparse
import warnings
import re
warnings.filterwarnings("ignore", category=SyntaxWarning)
from ete3 import NCBITaxa


def parse_list(list_of_taxids):
    """Parse a list of species/taxid and return dict {taxid: [species_name, ...]}"""
    result = {}
    for item in list_of_taxids:
        if '/' in item:
            species, taxid = item.rsplit('/', 1)
            taxid = int(taxid)
            if taxid in result:
                if isinstance(result[taxid], list):
                    result[taxid].append(species)
                else:
                    result[taxid] = [result[taxid], species]
            else:
                result[taxid] = species
    return result


def extract_taxid(parsed_list):
    """Return list of taxids"""
    return list(parsed_list.keys())


def download_newick(taxids):
    """Build a Newick tree from a list of taxids and rename leaves with Name|taxID"""
    ncbi = NCBITaxa()
    tree = ncbi.get_topology(taxids)
    return tree.write(format=1, format_root_node=True)

def rename_leafs(tree, name_map):
    """
    Remplace chaque taxid dans la chaîne Newick par le label Name|taxID.
    Si plusieurs noms pour un taxid, utilise le premier.
    """
    for taxid, names in name_map.items():
        if isinstance(names, list):
            label = f"{names[0]}|{taxid}"
        else:
            label = f"{names}|{taxid}"
        # Remplace le taxid suivi de ')' ou ',' ou ':' (feuille ou interne)
        tree = tree.replace(f"{taxid}", label)
    return tree


def find_nodes(tree):
    """
    Affiche tous les taxids (nœuds) trouvés dans la chaîne Newick,
    encadrés à gauche par '(' ou ',' et à droite par ',' ou ')' ou ';'
    """
    pattern = r'\)(\d+):'
    matches = re.findall(pattern, tree)
    for taxid in matches:
        name = get_name(taxid)
        label = f"{name}|{taxid}"
        tree = tree.replace(f"{taxid}", label)

    return tree


def get_name(taxid):
    ncbi = NCBITaxa()
    taxid_int = int(taxid)
    name = ncbi.get_taxid_translator([taxid_int])
    return name[taxid_int]

def main():
    parser = argparse.ArgumentParser(description="Build a Newick tree from a list of species/taxid")
    parser.add_argument("taxa", nargs="+", help="List of species/taxid, e.g. Homo_sapiens/9606")
    parser.add_argument("--with-distances", action="store_true", help="Add dummy branch lengths (=1)")
    parser.add_argument("--outfile", type=str, help="Output file (default: stdout)")

    args = parser.parse_args()

    name_map = parse_list(args.taxa)
    taxids = extract_taxid(name_map)
    newick = download_newick(taxids)
    newick = rename_leafs(newick, name_map)
    newick = find_nodes(newick)

    if not args.with_distances:
        newick = newick.replace(":1", "")

    if args.outfile:
        with open(args.outfile, "w") as f:
            f.write(newick + "\n")
        print(f"Newick tree saved to {args.outfile}")
    else:
        print(newick)

if __name__ == "__main__":
    main()
