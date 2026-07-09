#!/usr/bin/env python3

import argparse
import re
import csv
import gffutils
from collections import defaultdict


BAD_TERMS = re.compile(r"\b(hypothetical|unknown|putative)\b", re.IGNORECASE)
GO_RE = re.compile(r"GO:\d+")
IPR_RE = re.compile(r"IPR\d+")
PFAM_RE = re.compile(r"PF\d+")


def split_values(value, sep=","):
    if value in ("", "-", None):
        return []
    return [x.strip() for x in value.split(sep) if x.strip() and x.strip() != "-"]


def clean_description(stitle):
    return re.sub(r"\s+OS=.*$", "", stitle).strip()


def is_bad(description):
    return BAD_TERMS.search(description) is not None


def add_unique(attr, key, values):
    if not values:
        return

    if key not in attr:
        attr[key] = []

    old = list(attr[key])
    seen = set(old)

    for v in values:
        if v not in seen:
            old.append(v)
            seen.add(v)

    attr[key] = old


def parse_blast(blast_file):
    hits_by_query = {}

    with open(blast_file) as f:
        header = re.split(r"\s+", f.readline().strip())

        for line in f:
            if not line.strip():
                continue

            values = re.split(r"\s+", line.strip(), maxsplit=len(header) - 1)
            row = dict(zip(header, values))

            q = row["qseqid"]
            desc = clean_description(row["stitle"])

            hit = {
                "evalue": float(row["evalue"]),
                "bitscore": float(row["bitscore"]),
                "desc": desc,
                "sacc": row["sacc"],
                "bad": is_bad(desc),
            }

            hits_by_query.setdefault(q, []).append(hit)

    best_hits = {}

    for q, hits in hits_by_query.items():
        non_bad = [h for h in hits if not h["bad"]]
        candidates = non_bad if non_bad else hits
        best_hits[q] = sorted(candidates, key=lambda h: (h["evalue"], -h["bitscore"]))[0]

    return best_hits


def parse_interproscan(ipr_file):
    """
    InterProScan TSV columns:
    0 protein_id
    3 analysis
    4 signature_accession
    5 signature_description
    11 interpro_accession
    12 interpro_description
    13 GO terms
    14 pathways
    """
    annotations = defaultdict(lambda: defaultdict(set))

    with open(ipr_file) as f:
        for line in f:
            if not line.strip():
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 15:
                continue

            protein_id = cols[0]
            analysis = cols[3]
            sig_acc = cols[4]
            sig_desc = cols[5]
            ipr_acc = cols[11]
            ipr_desc = cols[12]
            go_col = cols[13]

            gos = GO_RE.findall(go_col)
            iprs = IPR_RE.findall(ipr_acc)
            pfams = PFAM_RE.findall(sig_acc)

            annotations[protein_id]["GO"].update(gos)
            annotations[protein_id]["InterPro"].update(iprs)
            annotations[protein_id]["PFAM"].update(pfams)

            if ipr_acc != "-" and ipr_desc != "-":
                annotations[protein_id]["Note"].add(f"{ipr_acc}:{ipr_desc}")

            if analysis in {"Phobius", "SignalP", "TMHMM"} and sig_acc != "-":
                annotations[protein_id]["Note"].add(f"{analysis}:{sig_acc}")

            if sig_desc != "-":
                annotations[protein_id]["domain"].add(sig_desc)

    return annotations


def parse_eggnog(eggnog_file):
    annotations = defaultdict(lambda: defaultdict(set))

    with open(eggnog_file) as f:
        header = None

        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue

            if not line.strip():
                continue

            cols = line.rstrip("\n").split("\t")
            row = dict(zip(header, cols))

            protein_id = row["query"]

            desc = row.get("Description", "-")
            name = row.get("Preferred_name", "-")

            gos = split_values(row.get("GOs", "-"), ",")
            ecs = split_values(row.get("EC", "-"), ",")
            kos = split_values(row.get("KEGG_ko", "-"), ",")
            pfams = split_values(row.get("PFAMs", "-"), ",")

            annotations[protein_id]["GO"].update(gos)
            annotations[protein_id]["EC_number"].update(ecs)
            annotations[protein_id]["KEGG_ko"].update(kos)
            annotations[protein_id]["PFAM"].update(pfams)

            if name != "-":
                annotations[protein_id]["Name"].add(name)

            if desc != "-" and not is_bad(desc):
                annotations[protein_id]["Note"].add(desc)

    return annotations


def annotate_gff(gff_file, output_file, blast_hits, ipr_annots, eggnog_annots):
    db = gffutils.create_db(
        gff_file,
        ":memory:",
        merge_strategy="create_unique",
        keep_order=True
    )

    with open(output_file, "w") as out:
        out.write("##gff-version 3\n")

        for feature in db.all_features(order_by=("seqid", "start")):

            if feature.featuretype == "mRNA":
                mrna_id = feature.id

                # BLASTp
                if mrna_id in blast_hits:
                    hit = blast_hits[mrna_id]
                    feature.attributes["product"] = [hit["desc"]]
                    add_unique(feature.attributes, "Dbxref", [f"UniProtKB:{hit['sacc']}"])
                    feature.attributes["blastp_evalue"] = [str(hit["evalue"])]
                    feature.attributes["blastp_bitscore"] = [str(hit["bitscore"])]

                # InterProScan
                if mrna_id in ipr_annots:
                    ann = ipr_annots[mrna_id]
                    add_unique(feature.attributes, "Ontology_term", sorted(ann["GO"]))
                    add_unique(feature.attributes, "Dbxref", [f"InterPro:{x}" for x in sorted(ann["InterPro"])])
                    add_unique(feature.attributes, "Dbxref", [f"PFAM:{x}" for x in sorted(ann["PFAM"])])
                    add_unique(feature.attributes, "Note", sorted(ann["Note"]))

                # eggNOG
                if mrna_id in eggnog_annots:
                    ann = eggnog_annots[mrna_id]
                    add_unique(feature.attributes, "Ontology_term", sorted(ann["GO"]))
                    add_unique(feature.attributes, "ec_number", sorted(ann["EC_number"]))
                    add_unique(feature.attributes, "Dbxref", sorted(ann["KEGG_ko"]))
                    add_unique(feature.attributes, "Dbxref", [f"PFAM:{x}" for x in sorted(ann["PFAM"])])
                    add_unique(feature.attributes, "Note", sorted(ann["Note"]))

                    if "Name" in ann and ann["Name"]:
                        feature.attributes["Name"] = [sorted(ann["Name"])[0]]

            out.write(str(feature) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Integrate BLASTp, InterProScan and eggNOG annotations into a GFF3 file."
    )

    parser.add_argument("--gff","-g", required=True)
    parser.add_argument("--output", "-o", required=True)
    parser.add_argument("--blast", "-b", required=False)
    parser.add_argument("--interproscan", "-i", required=False)
    parser.add_argument("--eggnog", "-e", required=False)

    args = parser.parse_args()

    blast_hits = parse_blast(args.blast) if args.blast else {}
    ipr_annots = parse_interproscan(args.interproscan) if args.interproscan else {}
    eggnog_annots = parse_eggnog(args.eggnog) if args.eggnog else {}

    annotate_gff(
        args.gff,
        args.output,
        blast_hits,
        ipr_annots,
        eggnog_annots
    )


if __name__ == "__main__":
    main()