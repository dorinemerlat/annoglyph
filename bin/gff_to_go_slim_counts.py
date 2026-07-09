#!/usr/bin/env python3

import argparse
import re
from collections import defaultdict, deque
from urllib.parse import unquote


GO_RE = re.compile(r"GO:\d+")


def parse_gff_attributes(attr_string):
    attrs = {}

    for field in attr_string.strip().split(";"):
        if not field:
            continue
        if "=" not in field:
            continue

        key, value = field.split("=", 1)
        attrs[key] = unquote(value)

    return attrs


def parse_obo(obo_file):
    """
    Parse GO OBO file.
    Returns:
      parents[go_id] = set(parent_go_ids)
      names[go_id] = term name
      obsolete = set(obsolete_go_ids)
    """
    parents = defaultdict(set)
    names = {}
    obsolete = set()

    current_id = None
    in_term = False

    with open(obo_file) as f:
        for line in f:
            line = line.rstrip("\n")

            if line == "[Term]":
                in_term = True
                current_id = None
                continue

            if line.startswith("[") and line != "[Term]":
                in_term = False
                current_id = None
                continue

            if not in_term:
                continue

            if line.startswith("id: "):
                current_id = line.replace("id: ", "").strip()
                continue

            if current_id is None:
                continue

            if line.startswith("name: "):
                names[current_id] = line.replace("name: ", "").strip()

            elif line.startswith("is_a: "):
                parent = line.split()[1]
                parents[current_id].add(parent)

            elif line.startswith("relationship: part_of "):
                parent = line.split()[2]
                parents[current_id].add(parent)

            elif line.startswith("is_obsolete: true"):
                obsolete.add(current_id)

    return parents, names, obsolete


def extract_goslim_terms(goslim_obo):
    goslim_terms = set()
    goslim_names = {}

    _, names, obsolete = parse_obo(goslim_obo)

    for go_id, name in names.items():
        if go_id not in obsolete:
            goslim_terms.add(go_id)
            goslim_names[go_id] = name

    return goslim_terms, goslim_names


def get_ancestors(go_id, parents):
    ancestors = set()
    queue = deque([go_id])

    while queue:
        current = queue.popleft()

        for parent in parents.get(current, []):
            if parent not in ancestors:
                ancestors.add(parent)
                queue.append(parent)

    return ancestors


def map_go_to_goslim(go_id, parents, goslim_terms):
    """
    A GO term maps to itself if it is already in GO slim,
    otherwise to any ancestor that belongs to GO slim.
    """
    mapped = set()

    if go_id in goslim_terms:
        mapped.add(go_id)

    ancestors = get_ancestors(go_id, parents)
    mapped.update(ancestors.intersection(goslim_terms))

    return mapped


def extract_go_per_mrna(gff_file):
    gene_to_go = defaultdict(set)

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue

            feature_type = cols[2]
            if feature_type != "mRNA":
                continue

            attrs = parse_gff_attributes(cols[8])

            mrna_id = attrs.get("ID")
            ontology = attrs.get("Ontology_term")

            if not mrna_id or not ontology:
                continue

            gos = GO_RE.findall(ontology)
            gene_to_go[mrna_id].update(gos)

    return gene_to_go


def write_goslim_counts(gene_to_go, parents, goslim_terms, goslim_names, output_file):
    goslim_to_genes = defaultdict(set)

    cache = {}

    for gene_id, go_terms in gene_to_go.items():
        gene_goslims = set()

        for go_id in go_terms:
            if go_id not in cache:
                cache[go_id] = map_go_to_goslim(go_id, parents, goslim_terms)

            gene_goslims.update(cache[go_id])

        for slim_id in gene_goslims:
            goslim_to_genes[slim_id].add(gene_id)

    with open(output_file, "w") as out:
        out.write("GO_slim\tname\tcount\n")

        for slim_id in sorted(
            goslim_to_genes,
            key=lambda x: (-len(goslim_to_genes[x]), x)
        ):
            name = goslim_names.get(slim_id, "NA")
            count = len(goslim_to_genes[slim_id])
            out.write(f"{slim_id}\t{name}\t{count}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate GO slim count table from a functionally annotated GFF3."
    )

    parser.add_argument("--gff", "-g", required=True, help="Annotated GFF3 file")
    parser.add_argument("--go-obo", required=True, help="go-basic.obo file")
    parser.add_argument("--goslim-obo", required=True, help="GO slim OBO file")
    parser.add_argument("--output", "-o", required=True, help="Output TSV file")

    args = parser.parse_args()

    parents, names, obsolete = parse_obo(args.go_obo)
    goslim_terms, goslim_names = extract_goslim_terms(args.goslim_obo)

    gene_to_go = extract_go_per_mrna(args.gff)

    write_goslim_counts(
        gene_to_go=gene_to_go,
        parents=parents,
        goslim_terms=goslim_terms,
        goslim_names=goslim_names,
        output_file=args.output
    )


if __name__ == "__main__":
    main()