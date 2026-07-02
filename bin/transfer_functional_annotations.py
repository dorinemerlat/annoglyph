#!/usr/bin/env python3
"""Transfer BLASTP, InterProScan and eggNOG-mapper annotations to a GFF3 file.

The script appends functional annotation attributes to mRNA/CDS features when IDs match
protein identifiers used by the external tools.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote


def parse_args():
    p = argparse.ArgumentParser(description="Add functional annotations to GFF3")
    p.add_argument("--gff", required=True)
    p.add_argument("--blast", required=True)
    p.add_argument("--interpro", required=True)
    p.add_argument("--eggnog", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def parse_attrs(s):
    d = {}
    for x in s.strip().split(";"):
        if not x:
            continue
        if "=" in x:
            k, v = x.split("=", 1)
            d[k] = v
    return d


def format_attrs(d):
    return ";".join(f"{k}={v}" for k, v in d.items() if v not in (None, ""))


def safe(v):
    return quote(str(v).replace(";", ","), safe="/:|,._- ")


def load_blast(path: Path):
    ann = {}
    if not path.exists():
        return ann
    with path.open(errors="ignore") as handle:
        header = None
        for line in handle:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if cols[0] == "qseqid":
                header = cols
                continue
            qid = cols[0]
            if qid in ann:
                continue
            if header and len(cols) == len(header):
                row = dict(zip(header, cols))
                ann[qid] = {
                    "blastp_hit": row.get("sseqid", ""),
                    "blastp_accession": row.get("sacc", row.get("saccver", "")),
                    "blastp_evalue": row.get("evalue", ""),
                    "blastp_bitscore": row.get("bitscore", row.get("bitscor", "")),
                    "blastp_description": row.get("stitle", ""),
                }
            elif len(cols) >= 12:
                ann[qid] = {
                    "blastp_hit": cols[1],
                    "blastp_evalue": cols[10],
                    "blastp_bitscore": cols[11],
                    "blastp_description": cols[-1] if len(cols) > 12 else "",
                }
    return ann


def load_interpro(path: Path):
    ann = defaultdict(lambda: defaultdict(set))
    if not path.exists():
        return ann
    with path.open(errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            qid = cols[0]
            analysis = cols[3]
            signature = cols[4]
            desc = cols[5] if len(cols) > 5 else ""
            ipr = cols[11] if len(cols) > 11 else ""
            ipr_desc = cols[12] if len(cols) > 12 else ""
            if signature and signature != "-":
                ann[qid]["interpro_signatures"].add(f"{analysis}:{signature}")
            if desc and desc != "-":
                ann[qid]["interpro_descriptions"].add(desc)
            if ipr and ipr != "-":
                ann[qid]["interpro_ids"].add(ipr)
            if ipr_desc and ipr_desc != "-":
                ann[qid]["interpro_descriptions"].add(ipr_desc)
    return ann


def load_eggnog(path: Path):
    ann = {}
    if not path.exists():
        return ann
    with path.open(errors="ignore") as handle:
        header = None
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#query") or line.startswith("query"):
                header = [x.lstrip("#") for x in line.rstrip("\n").split("\t")]
                continue
            cols = line.rstrip("\n").split("\t")
            if header and len(cols) == len(header):
                row = dict(zip(header, cols))
            else:
                row = {}
                if len(cols) >= 1:
                    row["query"] = cols[0]
                if len(cols) >= 8:
                    row["Description"] = cols[7]
                if len(cols) >= 10:
                    row["GOs"] = cols[9]
                if len(cols) >= 12:
                    row["KEGG_ko"] = cols[11]
            qid = row.get("query") or row.get("query_name") or cols[0]
            ann[qid] = {
                "eggnog_description": row.get("Description", row.get("description", "")),
                "eggnog_go": row.get("GOs", row.get("GO_terms", "")),
                "eggnog_ko": row.get("KEGG_ko", row.get("KEGG KO", "")),
                "eggnog_cog_category": row.get("COG_category", ""),
                "eggnog_preferred_name": row.get("Preferred_name", ""),
            }
    return ann


def candidate_ids(attrs):
    ids = []
    for key in ("ID", "Name", "protein_id", "transcript_id", "Parent"):
        if key in attrs:
            ids.extend(attrs[key].split(","))
    # Also try common transcript/protein suffix normalization.
    out = []
    for x in ids:
        out.append(x)
        out.append(x.replace(".t1", ".p1"))
        out.append(x.replace(".mRNA", ".protein"))
    return out


def main():
    args = parse_args()
    blast = load_blast(Path(args.blast))
    interpro = load_interpro(Path(args.interpro))
    eggnog = load_eggnog(Path(args.eggnog))

    with Path(args.gff).open(errors="ignore") as inp, Path(args.output).open("w") as out:
        for line in inp:
            if not line.strip() or line.startswith("#"):
                out.write(line)
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                out.write(line)
                continue
            attrs = parse_attrs(cols[8])
            keys = candidate_ids(attrs)
            merged = {}
            for k in keys:
                if k in blast:
                    merged.update(blast[k])
                if k in interpro:
                    for akey, vals in interpro[k].items():
                        merged[akey] = ",".join(sorted(vals))
                if k in eggnog:
                    merged.update(eggnog[k])
            for k, v in merged.items():
                if v and v != "-":
                    attrs[k] = safe(v)
            cols[8] = format_attrs(attrs)
            out.write("\t".join(cols) + "\n")


if __name__ == "__main__":
    main()
