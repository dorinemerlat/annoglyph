#!/usr/bin/env python3
"""
Generate TSV summary tables for annotation comparisons.

Outputs:
  annotation_summary.tsv
  protein_lengths.tsv

At least one input type must be provided:
  --proteins
  --busco
  --omark
  --gff

Filename format expected:
  ${genome_id}_${method_annotation}_${file_type}.${extension}
"""

from __future__ import annotations

import argparse
import json
import re
from functools import reduce
from pathlib import Path

import pandas as pd
from Bio import SeqIO
import gffutils


# -------------------------
# Path and filename parsing
# -------------------------

def collect_files(inputs, extensions):
    """
    Collect files from a list of paths or directories, filtering by extensions.
    """
    files = []

    for item in inputs:
        path = Path(item)

        if path.is_file():
            files.append(path)
        
        elif path.is_dir():
            for ext in extensions:
                files.extend(path.glob(f"*{ext}"))

        else:
            print(f"WARNING: input does not exist or matches nothing: {item}")

    files = sorted(set(files))
    return files



def parse_named_file(path):
    """
    Parse genome_id and method_annotation from:
      genome_id_method_annotation_filetype.ext

    Assumption:
      genome_id does not contain underscores.
      method_annotation may contain underscores.
    """
    
    file_name = path.stem

    genome_id= file_name.split("_")[0]
    method_annotation = file_name.split("_")[1]

    return {
        "genome_id": genome_id,
        "method_annotation": method_annotation,
    }


# -------------------------
# Protein parsing
# -------------------------

def parse_proteins(files):
    gene_counts = []
    protein_lengths = []

    for file in files:
        parsed = parse_named_file(file)
        genome_id = parsed["genome_id"]
        method = parsed["method_annotation"]

        print(f"[proteins] {genome_id} / {method}")

        n_proteins = 0
        lengths = []

        for record in SeqIO.parse(file, "fasta"):
            seq = str(record.seq).rstrip("*")
            protein_length = len(seq)

            protein_lengths.append({
                "genome_id": genome_id,
                "method_annotation": method,
                "protein_length": protein_length,
            })

            lengths.append(protein_length)
            n_proteins += 1

        mean_protein_length = sum(lengths) / len(lengths) if lengths else None
        q1_protein_length = pd.Series(lengths).quantile(0.25) if lengths else None

        gene_counts.append({
            "genome_id": genome_id,
            "method_annotation": method,
            "gene_number": n_proteins,
            "mean_protein_length": mean_protein_length,
            "q1_protein_length": q1_protein_length,
        })

    return pd.DataFrame(gene_counts), pd.DataFrame(protein_lengths)


# -------------------------
# BUSCO parsing
# -------------------------

def parse_busco(files):
    rows = []

    for file in files:
        parsed = parse_named_file(file)
        genome_id = parsed["genome_id"]
        method = parsed["method_annotation"]

        print(f"[BUSCO] {genome_id} / {method}")

        with open(file) as handle:
            data = json.load(handle)

        results = data.get("results", {})

        rows.append({
            "genome_id": genome_id,
            "method_annotation": method,
            "busco_complete_percent": results.get("Complete percentage"),
            "busco_single_percent": results.get("Single copy percentage"),
            "busco_duplicated_percent": results.get("Multi copy percentage"),
            "busco_fragmented_percent": results.get("Fragmented percentage"),
            "busco_missing_percent": results.get("Missing percentage"),
            "busco_complete_count": results.get("Complete BUSCOs"),
            "busco_single_count": results.get("Single copy BUSCOs"),
            "busco_duplicated_count": results.get("Multi copy BUSCOs"),
            "busco_fragmented_count": results.get("Fragmented BUSCOs"),
            "busco_missing_count": results.get("Missing BUSCOs"),
            "busco_n_markers": results.get("n_markers"),
            "busco_lineage": data.get("lineage_dataset", {}).get("name"),
        })

    return pd.DataFrame(rows)


# -------------------------
# OMArk parsing
# -------------------------

def parse_omark(files):
    rows = []

    for file in files:
        parsed = parse_named_file(file)

        row = {
            "genome_id": parsed["genome_id"],
            "method_annotation": parsed["method_annotation"],
        }

        print(f"[OMArk] {row['genome_id']} / {row['method_annotation']}")

        with open(file) as handle:
            for line in handle:
                line = line.strip()

                protein_match = re.search(
                    r"On the whole proteome, there \w+ ([0-9]+) proteins",
                    line,
                )
                if protein_match:
                    row["omark_protein_number"] = int(protein_match.group(1))

                hog_match = re.search(
                    r"S:([-0-9.]+)%,D:([-0-9.]+)%\[U:([-0-9.]+)%,E:([-0-9.]+)%\],M:([-0-9.]+)%",
                    line,
                )
                if hog_match:
                    row["omark_single_percent"] = float(hog_match.group(1))
                    row["omark_duplicated_percent"] = float(hog_match.group(2))
                    row["omark_unexpected_duplicated_percent"] = float(hog_match.group(3))
                    row["omark_expected_duplicated_percent"] = float(hog_match.group(4))
                    row["omark_missing_percent"] = float(hog_match.group(5))
                    row["omark_complete_percent"] = (
                        row["omark_single_percent"]
                        + row["omark_duplicated_percent"]
                    )

                consistency_match = re.search(
                    r"A:([-0-9.]+)%\[P:([-0-9.]+)%,F:([-0-9.]+)%\],"
                    r"I:([-0-9.]+)%\[P:([-0-9.]+)%,F:([-0-9.]+)%\],"
                    r"C:([-0-9.]+)%\[P:([-0-9.]+)%,F:([-0-9.]+)%\],"
                    r"U:([-0-9.]+)%",
                    line,
                )
                if consistency_match:
                    vals = [float(consistency_match.group(i)) for i in range(1, 11)]
                    a, ap, af, inc, ip, inf, c, cp, cf, u = vals

                    row.update({
                        "omark_consistent_percent": a,
                        "omark_consistent_partial_percent": ap,
                        "omark_consistent_fragmented_percent": af,
                        "omark_consistent_complete_percent": a - ap - af,

                        "omark_inconsistent_percent": inc,
                        "omark_inconsistent_partial_percent": ip,
                        "omark_inconsistent_fragmented_percent": inf,
                        "omark_inconsistent_complete_percent": inc - ip - inf,

                        "omark_contamination_percent": c,
                        "omark_contamination_partial_percent": cp,
                        "omark_contamination_fragmented_percent": cf,
                        "omark_contamination_complete_percent": c - cp - cf,

                        "omark_unknown_percent": u,
                    })

        rows.append(row)

    return pd.DataFrame(rows)


# -------------------------
# GFF parsing
# -------------------------

def get_attr(feature, key):
    value = feature.attributes.get(key, [])
    if isinstance(value, list):
        return value[0] if value else None
    return value


def has_product(feature):
    product = get_attr(feature, "product")
    return product is not None and str(product).strip() not in ["", ".", "NA"]


def has_go(feature):
    go = feature.attributes.get("Ontology_term", [])
    if isinstance(go, str):
        go = [go]
    return any(str(x).startswith("GO:") for x in go)


def get_canonical_transcripts(db):
    """
    Select one transcript per gene:
      1. longest cumulative CDS length
      2. longest cumulative exon length
      3. lexicographic transcript ID
    """

    canonical = {}

    for gene in db.features_of_type("gene"):
        transcripts = list(db.children(gene, featuretype="mRNA", order_by="start"))

        if not transcripts:
            continue

        scores = []

        for mrna in transcripts:
            cds_len = sum(
                len(cds)
                for cds in db.children(mrna, featuretype="CDS", order_by="start")
            )

            exon_len = sum(
                len(exon)
                for exon in db.children(mrna, featuretype="exon", order_by="start")
            )

            scores.append((mrna.id, cds_len, exon_len))

        best_id, _, _ = sorted(
            scores,
            key=lambda x: (-x[1], -x[2], x[0])
        )[0]

        canonical[gene.id] = best_id

    return canonical


def feature_length(feature):
    return abs(int(feature.end) - int(feature.start)) + 1


def parse_gff(files):
    if gffutils is None:
        raise ImportError(
            "gffutils is required for GFF parsing. Install with: pip install gffutils"
        )

    rows = []
    intron_lengths = []

    for gff in files:
        parsed = parse_named_file(gff)
        genome_id = parsed["genome_id"]
        method = parsed["method_annotation"]

        print(f"[GFF] {genome_id} / {method}")

        db = gffutils.create_db(
            str(gff),
            ":memory:",
            merge_strategy="create_unique",
            keep_order=True,
            force=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True,
        )

        canonical_ids = set(get_canonical_transcripts(db).values())
        n_genes = len(canonical_ids)

        n_start_and_stop = 0
        n_start_only = 0
        n_stop_only = 0

        exon_counts = []
        n_monoexonic = 0

        n_product_only = 0
        n_go_only = 0
        n_product_and_go = 0

        seen_introns = set()

        for mrna_id in canonical_ids:
            mrna = db[mrna_id]

            has_start = bool(list(
                db.children(mrna, featuretype="start_codon", order_by="start")
            ))

            has_stop = bool(list(
                db.children(mrna, featuretype="stop_codon", order_by="start")
            ))

            if has_start and has_stop:
                n_start_and_stop += 1
            elif has_start:
                n_start_only += 1
            elif has_stop:
                n_stop_only += 1

            exons = list(db.children(mrna, featuretype="exon", order_by="start"))
            n_exons = len(exons)

            exon_counts.append(n_exons)
            n_monoexonic += int(n_exons == 1)

            introns = list(db.children(mrna, featuretype="intron", order_by="start"))

            for intron in introns:
                intron_id = str(intron.id)

                # Ignore AGAT-added duplicated introns if present
                if intron_id.startswith("intron_added"):
                    continue

                intron_key = (
                    mrna_id,
                    intron.seqid,
                    int(intron.start),
                    int(intron.end),
                    intron.strand,
                )

                if intron_key in seen_introns:
                    continue

                seen_introns.add(intron_key)

                intron_lengths.append({
                    "genome_id": genome_id,
                    "method_annotation": method,
                    "transcript_id": mrna_id,
                    "intron_length": feature_length(intron),
                })

            product = has_product(mrna)
            go = has_go(mrna)

            if product and go:
                n_product_and_go += 1
            elif product:
                n_product_only += 1
            elif go:
                n_go_only += 1

        if n_genes == 0:
            rows.append({
                "genome_id": genome_id,
                "method_annotation": method,
                "gff_gene_number": 0,
                "mean_exons_per_gene": None,
                "percent_monoexonic_genes": None,
                "percent_genes_start_and_stop": None,
                "percent_genes_start_only": None,
                "percent_genes_stop_only": None,
                "percent_genes_product_only": None,
                "percent_genes_go_only": None,
                "percent_genes_product_and_go": None,
            })
            continue

        rows.append({
            "genome_id": genome_id,
            "method_annotation": method,
            "gff_gene_number": n_genes,

            "percent_genes_start_and_stop": n_start_and_stop / n_genes * 100,
            "percent_genes_start_only": n_start_only / n_genes * 100,
            "percent_genes_stop_only": n_stop_only / n_genes * 100,

            "mean_exons_per_gene": sum(exon_counts) / n_genes,
            "percent_monoexonic_genes": n_monoexonic / n_genes * 100,

            "percent_genes_product_only": n_product_only / n_genes * 100,
            "percent_genes_go_only": n_go_only / n_genes * 100,
            "percent_genes_product_and_go": n_product_and_go / n_genes * 100,
        })

    return pd.DataFrame(rows), pd.DataFrame(intron_lengths, columns=["genome_id", "method_annotation", "transcript_id", "intron_length"])


# -------------------------
# Summary merging
# -------------------------

def merge_summary_tables(tables):
    tables = [
        df for df in tables
        if not df.empty and {"genome_id", "method_annotation"}.issubset(df.columns)
    ]

    if not tables:
        return pd.DataFrame()

    summary = reduce(
        lambda left, right: pd.merge(
            left,
            right,
            on=["genome_id", "method_annotation"],
            how="outer",
        ),
        tables,
    )

    return summary.sort_values(["genome_id", "method_annotation"])


def save_tsv(df, path):
    df.to_csv(path, sep="\t", index=False, float_format="%.6g")
    print(f"Saved {path}")


# -------------------------
# Main
# -------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate annotation comparison TSV tables."
    )

    parser.add_argument(
        "-p",
        "--proteins",
        nargs="+",
        default=None,
        help="Protein FASTA files, directories, or glob patterns.",
    )

    parser.add_argument(
        "-b",
        "--busco",
        nargs="+",
        default=None,
        help="BUSCO JSON files, directories, or glob patterns.",
    )

    parser.add_argument(
        "-m",
        "--omark",
        nargs="+",
        default=None,
        help="OMArk .sum files, directories, or glob patterns.",
    )

    parser.add_argument(
        "-g",
        "--gff",
        nargs="+",
        default=None,
        help="GFF/GFF3 files, directories, or glob patterns.",
    )

    parser.add_argument(
        "-o",
        "--out",
        required=True,
        type=Path,
        help="Base name for output files",
    )

    parser.add_argument("--genome-id")
    parser.add_argument("--annotation-method")
    args = parser.parse_args()

    if not any([args.proteins, args.busco, args.omark, args.gff]):
        parser.error("At least one input argument is required: --proteins, --busco, --omark or --gff")


    protein_files = collect_files(args.proteins, (".fasta", ".fa", ".faa"))
    busco_files = collect_files(args.busco, (".json",))
    omark_files = collect_files(args.omark, (".sum",))
    gff_files = collect_files(args.gff, (".gff", ".gff3"))

    df_gene_number = pd.DataFrame()
    df_protein_lengths = pd.DataFrame()
    df_busco = pd.DataFrame()
    df_omark = pd.DataFrame()
    df_gff = pd.DataFrame()
    df_intron_lengths = pd.DataFrame()

    if protein_files:
        df_gene_number, df_protein_lengths = parse_proteins(protein_files)

    if busco_files:
        df_busco = parse_busco(busco_files)

    if omark_files:
        df_omark = parse_omark(omark_files)

    if gff_files:
        df_gff, df_intron_lengths = parse_gff(gff_files)

    # Explicit metadata avoids relying on underscores in staged filenames.
    for table in [df_gene_number, df_protein_lengths, df_busco, df_omark, df_gff, df_intron_lengths]:
        if not table.empty:
            if args.genome_id:
                table["genome_id"] = args.genome_id
            if args.annotation_method:
                table["method_annotation"] = args.annotation_method

    summary = merge_summary_tables([
        df_gene_number,
        df_busco,
        df_omark,
        df_gff,
    ])

    save_tsv(summary, Path(str(args.out) + "_annotation_summary.tsv"))

    save_tsv(df_protein_lengths, Path(str(args.out) + "_protein_lengths.tsv"))

    save_tsv(df_intron_lengths, Path(str(args.out) + "_intron_lengths.tsv"))

if __name__ == "__main__":
    main()