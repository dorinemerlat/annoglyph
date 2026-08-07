#!/usr/bin/env python3

import argparse
import re
import sys
from urllib.parse import unquote

import gffutils
from Bio import SeqIO
from Bio.Seq import Seq

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Transfer functional annotations from GFF3 mRNA features "
            "to protein FASTA headers."
        )
    )
    parser.add_argument("-f", "--fasta", required=True, help="Input protein FASTA file.")
    parser.add_argument("-g", "--gff", required=True, help="Input GFF3 annotation file.")
    parser.add_argument("-o", "--output", required=True, help="Output annotated protein FASTA file.")
    parser.add_argument("-s", "--organism", help="Organism name to add as OS= in FASTA headers.")
    parser.add_argument("-t", "--taxid", help="NCBI taxonomy identifier to add as OX= in FASTA headers.")
    return parser.parse_args()


def first_attribute(feature, attribute_names):
    """
    Return the first non-empty value found among the requested
    GFF3 attributes.
    """

    for attribute_name in attribute_names:
        values = feature.attributes.get(attribute_name, [])

        if not values:
            continue

        value = unquote(str(values[0])).strip()

        if value:
            return value

    return None


def clean_text(text):
    """
    Clean text before inserting it into a one-line FASTA header.
    """

    if not text:
        return None

    text = unquote(text)
    text = re.sub(r"\s+", " ", text).strip()
    text = text.rstrip(";,")
    return text or None


def get_product(feature):
    """
    Determine the protein name using the following priority:

        product
        Name
        Note
        Uncharacterized protein
    """

    product = first_attribute(feature, ["product", "protein_name"])

    if product:
        return clean_text(product)

    name = first_attribute(feature, ["Name"])

    if name:
        return clean_text(name)

    note = first_attribute(feature, ["Note"])

    if note:
        note = clean_text(note)

        # Notes often start with an InterPro accession:
        # IPR026983:Dynein heavy chain,activity...
        note = re.sub(r"^(?:IPR\d+):", "", note)

        return note

    return "Uncharacterized protein"


def get_gene_id(feature):
    """
    Retrieve the parent gene ID from the mRNA Parent attribute.
    """

    gene_id = first_attribute(feature, ["gene", "gene_id", "Parent"])

    if gene_id:
        return clean_text(gene_id)

    return None


def build_annotation_index(db, organism, taxid):
    """
    Build a dictionary indexed by mRNA ID.

    The organism name (OS) and NCBI TaxID (OX) are provided by the
    pipeline and therefore are not read from the GFF file.

    Example:
        annotations["STRAC3.g4.t1"] = {
            "product": "Dynein heavy chain, cytoplasmic",
            "gene_id": "STRAC3.g4",
            "organism": "Strigamia acuminata",
            "taxid": "1255758"
        }
    """

    organism = clean_text(organism)
    taxid = clean_text(str(taxid))

    if not organism:
        raise ValueError("No organism name was provided.")

    if not taxid:
        raise ValueError("No NCBI TaxID was provided.")

    annotations = {}

    for feature in db.features_of_type("mRNA", order_by=("seqid", "start")):
        transcript_id = feature.id

        if not transcript_id:
            transcript_id = first_attribute(
                feature,
                ["ID", "transcript_id"],
            )

        if not transcript_id:
            continue

        transcript_id = clean_text(transcript_id)

        annotations[transcript_id] = {
            "product": get_product(feature),
            "gene_id": get_gene_id(feature),
            "organism": organism,
            "taxid": taxid,
        }

    return annotations


def candidate_fasta_ids(record):
    """
    Generate possible IDs for matching a FASTA sequence against
    a GFF mRNA ID.
    """

    candidates = [
        record.id,
        record.name,
        record.description.split()[0],
    ]

    expanded = []

    for candidate in candidates:
        candidate = candidate.strip()

        if not candidate:
            continue

        expanded.append(candidate)

        # Handle possible prefixes.
        expanded.append(
            re.sub(
                r"^(?:transcript:|mRNA:|protein:)",
                "",
                candidate,
            )
        )

        # Handle pipe-delimited headers.
        if "|" in candidate:
            expanded.extend(candidate.split("|"))

    # Remove duplicates while preserving order.
    return list(dict.fromkeys(expanded))


def find_annotation(record, annotations):
    """
    Find the GFF annotation corresponding to a FASTA record.
    """

    for candidate in candidate_fasta_ids(record):
        if candidate in annotations:
            return candidate, annotations[candidate]

    return None, None


def make_header(sequence_id, annotation):
    """
    Construct a UniProt-inspired FASTA header.

    Example:
        STRAC3.g4.t1 Dynein heavy chain, cytoplasmic GN=STRAC3.g4
    """

    parts = [
        sequence_id,
        annotation["product"],
    ]

    if annotation["organism"]:
        parts.append(f"OS={annotation['organism']}")

    if annotation["taxid"]:
        parts.append(f"OX={annotation['taxid']}")

    if annotation["gene_id"]:
        parts.append(f"GN={annotation['gene_id']}")

    return " ".join(parts)


def main():
    args = parse_args()

    print(f"Reading GFF: {args.gff}")

    try:
        db = gffutils.create_db(
            args.gff,
            ":memory:",
            merge_strategy="create_unique",
            keep_order=True,
        )
    except Exception as error:
        sys.exit(f"ERROR: could not create GFF database: {error}")

    annotations = build_annotation_index(db, args.organism, args.taxid)

    if not annotations:
        sys.exit(
            "ERROR: no mRNA annotations were found in the GFF. "
            "Check that the feature type is 'mRNA'."
        )

    print(f"mRNA annotations loaded: {len(annotations)}")
    print(f"Reading FASTA: {args.fasta}")

    n_proteins = 0
    n_annotated = 0
    n_unmatched = 0
    unmatched_ids = []

    try:
        input_records = SeqIO.parse(args.fasta, "fasta")

        with open(args.output, "w", encoding="utf-8") as output_handle:
            for record in input_records:
                n_proteins += 1

                matched_id, annotation = find_annotation(
                    record,
                    annotations,
                )

                if annotation is None:
                    sequence_id = record.id

                    fallback_annotation = {
                        "product": "Uncharacterized protein",
                        "gene_id": None,
                        "organism": clean_text(args.organism),
                        "taxid": clean_text(str(args.taxid)),
                    }

                    header = make_header(
                        sequence_id,
                        fallback_annotation,
                    )

                    n_unmatched += 1
                    unmatched_ids.append(record.id)

                else:
                    sequence_id = matched_id
                    header = make_header(
                        sequence_id,
                        annotation,
                    )

                    n_annotated += 1

                # When the description begins with record.id, Biopython
                # writes the provided description as the complete header
                # without duplicating the identifier.
                record.id = sequence_id
                record.name = sequence_id
                record.description = header
                record.seq = Seq(str(record.seq).rstrip("*"))

                SeqIO.write(
                    record,
                    output_handle,
                    "fasta",
                )

    except Exception as error:
        sys.exit(f"ERROR while processing FASTA: {error}")

    print()
    print(f"Proteins processed : {n_proteins}")
    print(f"Proteins annotated : {n_annotated}")
    print(f"Proteins unmatched : {n_unmatched}")
    print(f"Written            : {args.output}")

    if unmatched_ids:
        unmatched_file = f"{args.output}.unmatched_ids.txt"

        with open(unmatched_file, "w", encoding="utf-8") as handle:
            for sequence_id in unmatched_ids:
                handle.write(f"{sequence_id}\n")

        print(f"Unmatched IDs      : {unmatched_file}")


if __name__ == "__main__":
    main()