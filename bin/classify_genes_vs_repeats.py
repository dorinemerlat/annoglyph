#!/usr/bin/env python3

import argparse
import subprocess
import tempfile
import gffutils
import sys
from collections import defaultdict


KEEP_REPEAT_TYPES = {"transposable_element", "unknown"}


def run(cmd, stdout=None):
    subprocess.run(cmd, check=True, stdout=stdout)


def sort_bed(path):
    subprocess.run(["sort", "-k1,1", "-k2,2n", path, "-o", path], check=True)


def parse_args():
    parser = argparse.ArgumentParser(description="Classify overlaps between protein-coding genes and repeats from GFF files.")
    parser.add_argument("--proteins", required=True, help="GFF with protein-coding gene annotations")
    parser.add_argument("--repeats", required=True, help="GFF with repeat annotations")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--border-size", type=int, default=100, help="Border size in bp around each repeat extremity (default: 100)")
    parser.add_argument("--merge-repeats", dest="merge_repeats", action="store_true", default=False, help="Merge repeats before classification (disabled by default)")
    parser.add_argument("--tmp", default=None, help="Temporary directory to use (optional)")
    return parser.parse_args()


def extract_genes(gff, bed_out):
    db = gffutils.create_db(
        gff,
        ":memory:",
        merge_strategy="create_unique",
        keep_order=True
    )

    with open(bed_out, "w") as out:
        for gene in db.features_of_type("gene"):
            gstart = int(gene.start)
            gend = int(gene.end)
            if gend < gstart:
                sys.stderr.write(f"Warning: in {gff} gene {gene.id} has end<{gstart} (end={gend}), swapping coordinates\n")
                gstart, gend = gend, gstart
            start0 = max(0, gstart - 1)
            out.write(f"{gene.chrom}\t{start0}\t{gend}\t{gene.id}\n")


def extract_repeats(gff, bed_out):
    """
    Output BED:
    chrom  start  end  repeat_id  repeat_label

    repeat_label:
      - if repeat_type == unknown:
            Unknown
      - if repeat_type == transposable_element:
            repeat_class/repeat_order
            and if one is NA -> Unknown
    """
    db = gffutils.create_db(
        gff,
        ":memory:",
        merge_strategy="create_unique",
        keep_order=True
    )

    with open(bed_out, "w") as out:
        for rep in db.all_features():
            attrs = rep.attributes

            repeat_type = attrs.get("repeat_type", ["NA"])[0]
            if repeat_type not in KEEP_REPEAT_TYPES:
                continue

            repeat_id = attrs.get("ID", [rep.id if rep.id else "."])[0]
            repeat_order = attrs.get("repeat_order", ["NA"])[0]
            repeat_class = attrs.get("repeat_class", ["NA"])[0]

            if repeat_type == "unknown":
                repeat_label = "Unknown"
            else:
                if repeat_class == "NA" or repeat_order == "NA":
                    repeat_label = "Unknown"
                else:
                    repeat_label = f"{repeat_class}/{repeat_order}"

            try:
                rstart = int(rep.start)
                rend = int(rep.end)
            except Exception:
                # fallback to original values if not parseable
                rstart = rep.start
                rend = rep.end

            if isinstance(rstart, int) and isinstance(rend, int) and rend < rstart:
                sys.stderr.write(f"Warning: in {gff} repeat {repeat_id} has end<{rstart} (end={rend}), swapping coordinates\n")
                rstart, rend = rend, rstart

            rstart0 = max(0, int(rstart) - 1)
            out.write(
                f"{rep.chrom}\t{rstart0}\t{rend}\t{repeat_id}\t{repeat_label}\n"
            )


def merge_repeats(repeats_bed, merged_bed):
    with open(merged_bed, "w") as out:
        run(["bedtools", "merge", "-i", repeats_bed], stdout=out)


def copy_repeat_coords(repeats_bed, out_bed):
    with open(repeats_bed) as fin, open(out_bed, "w") as out:
        for line in fin:
            p = line.rstrip("\n").split("\t")
            out.write(f"{p[0]}\t{p[1]}\t{p[2]}\n")


def read_gene_ids(genes_bed):
    ids = set()
    with open(genes_bed) as f:
        for line in f:
            ids.add(line.rstrip("\n").split("\t")[3])
    return ids


def intersect_lines(a_bed, b_bed):
    res = subprocess.check_output(
        ["bedtools", "intersect", "-a", a_bed, "-b", b_bed, "-wa", "-wb"],
        text=True
    )
    return [line for line in res.splitlines() if line.strip()]


def get_gene_inside_repeat(genes_bed, repeats_coords_bed):
    """
    gene fully inside repeat
    genes:   4 cols
    repeats: 3 cols
    """
    gene_to_repeat_coords = defaultdict(set)

    for line in intersect_lines(genes_bed, repeats_coords_bed):
        p = line.split("\t")
        gene_id = p[3]
        gstart = int(p[1])
        gend = int(p[2])
        rchrom = p[4]
        rstart = int(p[5])
        rend = int(p[6])

        if gstart >= rstart and gend <= rend:
            gene_to_repeat_coords[gene_id].add((rchrom, rstart, rend))

    return gene_to_repeat_coords


def get_repeat_inside_gene(genes_bed, repeats_coords_bed):
    """
    repeat fully inside gene
    """
    gene_to_repeat_coords = defaultdict(set)

    for line in intersect_lines(genes_bed, repeats_coords_bed):
        p = line.split("\t")
        gene_id = p[3]
        gstart = int(p[1])
        gend = int(p[2])
        rchrom = p[4]
        rstart = int(p[5])
        rend = int(p[6])

        if rstart >= gstart and rend <= gend:
            gene_to_repeat_coords[gene_id].add((rchrom, rstart, rend))

    return gene_to_repeat_coords


def get_overlap_states(genes_bed, repeats_coords_bed, border_size):
    """
    Returns for each gene:
      - repeat overlaps
      - border overlaps

    Border definition:
      for each repeat [rstart, rend):
        left border  = [max(0, rstart-border_size), rstart+border_size)
        right border = [max(0, rend-border_size), rend+border_size)

    Categories are assigned later:
      internal_overlap: overlap repeat AND overlap border
      external_overlap: overlap border AND NOT overlap repeat
      other_overlap: overlap repeat AND NOT overlap border
    """
    gene_to_states = defaultdict(lambda: {
        "repeat": set(),
        "border": set()
    })

    for line in intersect_lines(genes_bed, repeats_coords_bed):
        p = line.split("\t")

        gene_id = p[3]
        gstart = int(p[1])
        gend = int(p[2])
        rchrom = p[4]
        rstart = int(p[5])
        rend = int(p[6])

        coord = (rchrom, rstart, rend)
        gene_to_states[gene_id]["repeat"].add(coord)

        left_start = max(0, rstart - border_size)
        left_end = rstart + border_size
        right_start = max(0, rend - border_size)
        right_end = rend + border_size

        overlaps_left_border = (gstart < left_end and gend > left_start)
        overlaps_right_border = (gstart < right_end and gend > right_start)

        if overlaps_left_border or overlaps_right_border:
            gene_to_states[gene_id]["border"].add(coord)

    return gene_to_states


def build_repeat_info_no_merge(genes_bed, repeats_bed):
    """
    Map for no-merge mode:
      gene_id -> {(chrom, start, end): set((repeat_id, repeat_label))}
    repeats_bed columns:
      chrom start end repeat_id repeat_label
    """
    info = defaultdict(lambda: defaultdict(set))

    for line in intersect_lines(genes_bed, repeats_bed):
        p = line.split("\t")
        gene_id = p[3]
        repeat_chrom = p[4]
        repeat_start = int(p[5])
        repeat_end = int(p[6])
        repeat_id = p[7]
        repeat_label = p[8]

        info[gene_id][(repeat_chrom, repeat_start, repeat_end)].add((repeat_id, repeat_label))

    return info


def build_repeat_info_with_merge(genes_bed, repeats_bed, repeats_coords_bed):
    """
    Map for merge mode:
      gene_id -> {(merged_chrom, merged_start, merged_end): set((repeat_id, repeat_label))}

    This fixes the issue where merged coordinates no longer match exact original repeat coordinates.
    """
    info = defaultdict(lambda: defaultdict(set))

    # Intersect original repeats against merged blocks to know which original repeats belong to each merged block
    # A = original repeats (5 cols)
    # B = merged coords    (3 cols)
    res = subprocess.check_output(
        ["bedtools", "intersect", "-a", repeats_bed, "-b", repeats_coords_bed, "-wa", "-wb"],
        text=True
    )

    merged_to_original = defaultdict(set)

    for line in res.splitlines():
        if not line.strip():
            continue

        p = line.split("\t")
        orig_repeat_id = p[3]
        orig_repeat_label = p[4]
        merged_coord = (p[5], int(p[6]), int(p[7]))
        merged_to_original[merged_coord].add((orig_repeat_id, orig_repeat_label))

    # Now intersect genes with merged coords to attach merged block -> original repeat annotations for each gene
    res2 = subprocess.check_output(
        ["bedtools", "intersect", "-a", genes_bed, "-b", repeats_coords_bed, "-wa", "-wb"],
        text=True
    )

    for line in res2.splitlines():
        if not line.strip():
            continue

        p = line.split("\t")
        gene_id = p[3]
        merged_coord = (p[4], int(p[5]), int(p[6]))

        for rep_id, rep_label in merged_to_original.get(merged_coord, set()):
            info[gene_id][merged_coord].add((rep_id, rep_label))

    return info


def flatten_repeat_info(gene_id, repeat_coords, repeat_info):
    repeat_ids = set()
    repeat_labels = set()

    for coord in repeat_coords:
        for rep_id, rep_label in repeat_info.get(gene_id, {}).get(coord, set()):
            repeat_ids.add(rep_id)
            repeat_labels.add(rep_label)

    if not repeat_ids:
        return ".", "."

    return ",".join(sorted(repeat_ids)), ",".join(sorted(repeat_labels))


def main():
    args = parse_args()

    with tempfile.TemporaryDirectory(dir=args.tmp) as tmp:
        genes_bed = f"{tmp}/genes.bed"
        repeats_bed = f"{tmp}/repeats.bed"
        repeats_coords_bed = f"{tmp}/repeats.coords.bed"

        extract_genes(args.proteins, genes_bed)
        extract_repeats(args.repeats, repeats_bed)

        sort_bed(genes_bed)
        sort_bed(repeats_bed)

        if args.merge_repeats:
            merge_repeats(repeats_bed, repeats_coords_bed)
        else:
            copy_repeat_coords(repeats_bed, repeats_coords_bed)
            sort_bed(repeats_coords_bed)

        all_genes = read_gene_ids(genes_bed)

        gene_inside = get_gene_inside_repeat(genes_bed, repeats_coords_bed)
        repeat_inside = get_repeat_inside_gene(genes_bed, repeats_coords_bed)
        overlap_states = get_overlap_states(genes_bed, repeats_coords_bed, args.border_size)

        if args.merge_repeats:
            repeat_info = build_repeat_info_with_merge(genes_bed, repeats_bed, repeats_coords_bed)
        else:
            repeat_info = build_repeat_info_no_merge(genes_bed, repeats_bed)

        assigned = {}
        results = []

        # 1. gene_inside_repeat
        for gene_id in sorted(gene_inside):
            ids, labels = flatten_repeat_info(gene_id, gene_inside[gene_id], repeat_info)
            assigned[gene_id] = "gene_inside_repeat"
            results.append((gene_id, "gene_inside_repeat", ids, labels))

        # 2. repeat_inside_gene
        for gene_id in sorted(repeat_inside):
            if gene_id in assigned:
                continue
            ids, labels = flatten_repeat_info(gene_id, repeat_inside[gene_id], repeat_info)
            assigned[gene_id] = "repeat_inside_gene"
            results.append((gene_id, "repeat_inside_gene", ids, labels))

        # 3/4/5. overlap-based classes
        for gene_id in sorted(all_genes):
            if gene_id in assigned:
                continue

            repeat_coords = overlap_states.get(gene_id, {}).get("repeat", set())
            border_coords = overlap_states.get(gene_id, {}).get("border", set())

            if repeat_coords and border_coords:
                class_name = "internal_overlap"
                coords = repeat_coords & border_coords
                if not coords:
                    coords = repeat_coords
                ids, labels = flatten_repeat_info(gene_id, coords, repeat_info)
                assigned[gene_id] = class_name
                results.append((gene_id, class_name, ids, labels))

            elif (not repeat_coords) and border_coords:
                class_name = "external_overlap"
                coords = border_coords
                ids, labels = flatten_repeat_info(gene_id, coords, repeat_info)
                assigned[gene_id] = class_name
                results.append((gene_id, class_name, ids, labels))

            elif repeat_coords and (not border_coords):
                class_name = "other_overlap"
                coords = repeat_coords
                ids, labels = flatten_repeat_info(gene_id, coords, repeat_info)
                assigned[gene_id] = class_name
                results.append((gene_id, class_name, ids, labels))

        # 6. no_overlap
        for gene_id in sorted(all_genes):
            if gene_id not in assigned:
                results.append((gene_id, "no_overlap", ".", "."))

        with open(args.output, "w") as out:
            out.write("gene_id\tclass\trepeat_id\trepeat_order\n")
            for gene_id, class_name, rep_ids, rep_labels in sorted(results, key=lambda x: x[0]):
                out.write(f"{gene_id}\t{class_name}\t{rep_ids}\t{rep_labels}\n")


if __name__ == "__main__":
    main()