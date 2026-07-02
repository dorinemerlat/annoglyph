#!/usr/bin/env python3

import argparse
import gffutils
import os


def parse_args():
    parser = argparse.ArgumentParser(description="Reformat the BRAKER GFF3 file to make it more compatible with OMARK scripts")
    parser.add_argument("-i", "--input", required=True, help="Input GFF3 file")
    parser.add_argument("-o", "--output", required=True, help="Output reformatted GFF3 file")
    return parser.parse_args()

def add_prefix(identifier, prefix):
    if identifier.startswith(prefix):
        return identifier
    return f"{prefix}{identifier}"

def strip_prefix(value, prefixes=None):
    if prefixes is None:
        prefixes = ("gene:", "transcript:")
    for prefix in prefixes:
        if value.startswith(prefix):
            return value[len(prefix):]
    return value

def copy_attributes(attrs):
    return {key: list(values) for key, values in attrs.items()}

def main():
    args = parse_args()
    input = args.input
    output = args.output

    db = gffutils.create_db(input, ':memory:', merge_strategy="create_unique", keep_order=True)

    transcript_prefix = "transcript:"
    gene_prefix = "gene:"
    
    # Premier passage : construire les correspondances d'identifiants
    with open(output, "w") as out:
        out.write("##gff-version 3\n")
        for line in db.all_features(order_by=("seqid", "start")):
            attributes = copy_attributes(line.attributes)
            if line.featuretype == "gene":
                old_gene_id = line.id
                new_gene_id = add_prefix(old_gene_id, gene_prefix)
                attributes ["ID"] = [new_gene_id]

            elif line.featuretype == "mRNA":
                old_transcript_id = line.id
                new_transcript_id = add_prefix(old_transcript_id, transcript_prefix)
                attributes["ID"] = [new_transcript_id]

                old_parent = attributes["Parent"][0]
                new_parent = add_prefix(old_parent, gene_prefix)
                attributes["Parent"] = [new_parent]
            else:
                old_parent = attributes["Parent"][0]
                new_parent = add_prefix(old_parent, transcript_prefix)
                attributes["Parent"] = [new_parent]

                if line.featuretype == "CDS":
                    protein_id = old_parent
                    attributes["protein_id"] = [protein_id]
                    attributes["Name"] = [protein_id]

            line.attributes = attributes
            out.write(str(line) + "\n")

if __name__ == "__main__":
    main()