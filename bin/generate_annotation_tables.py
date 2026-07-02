#!/usr/bin/env python3
"""Generate compact annotation QC tables from GFF, protein FASTA, BUSCO JSON and OMARK summary.

Outputs three TSV files using the provided prefix:
  <prefix>_summary.tsv
  <prefix>_protein_lengths.tsv
  <prefix>_intron_lengths.tsv
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path
from statistics import mean, median
from typing import Dict, Iterable, List, Tuple


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate annotation comparison tables")
    p.add_argument("--gff", required=True, help="Input GFF/GFF3 file")
    p.add_argument("--proteins", "--fasta", dest="proteins", required=True, help="Canonical protein FASTA")
    p.add_argument("--busco", "--busco-json", dest="busco", required=True, help="BUSCO JSON summary")
    p.add_argument("--omark", "--omark-sum", dest="omark", required=True, help="OMArk .sum file")
    p.add_argument("--prefix", required=True, help="Output prefix")
    return p.parse_args()


def fasta_lengths(path: Path) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    name = None
    seq_len = 0
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = seq_len
                name = line[1:].split()[0]
                seq_len = 0
            else:
                seq_len += len(line)
    if name is not None:
        lengths[name] = seq_len
    return lengths


def parse_attributes(attr: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for field in attr.strip().split(";"):
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
            out[k] = v
        elif " " in field:
            k, v = field.split(" ", 1)
            out[k] = v.strip('"')
    return out


def parse_gff(path: Path) -> dict:
    genes = {}
    mrna_to_gene = {}
    exon_by_tx = defaultdict(list)
    cds_by_tx = defaultdict(list)
    features_count = defaultdict(int)

    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue
            attr = parse_attributes(attrs)
            features_count[ftype] += 1
            if ftype == "gene":
                gid = attr.get("ID") or attr.get("gene_id")
                if gid:
                    genes[gid] = (seqid, start_i, end_i, strand)
            elif ftype in {"mRNA", "transcript"}:
                tid = attr.get("ID") or attr.get("transcript_id")
                parent = attr.get("Parent") or attr.get("gene_id")
                if tid and parent:
                    mrna_to_gene[tid] = parent.split(",")[0]
            elif ftype == "exon":
                parent = attr.get("Parent") or attr.get("transcript_id")
                if parent:
                    for tid in parent.split(","):
                        exon_by_tx[tid].append((start_i, end_i))
            elif ftype == "CDS":
                parent = attr.get("Parent") or attr.get("transcript_id")
                if parent:
                    for tid in parent.split(","):
                        cds_by_tx[tid].append((start_i, end_i))

    exon_counts = [len(v) for v in exon_by_tx.values()]
    cds_counts = [len(v) for v in cds_by_tx.values()]
    intron_lengths = []
    for tx, exons in exon_by_tx.items():
        if len(exons) < 2:
            continue
        exons = sorted(exons)
        for (_, prev_end), (next_start, _) in zip(exons, exons[1:]):
            intron = next_start - prev_end - 1
            if intron > 0:
                intron_lengths.append(intron)

    mono_exonic = sum(1 for c in exon_counts if c == 1)
    n_transcripts = len(set(mrna_to_gene) | set(exon_by_tx) | set(cds_by_tx))
    n_genes = len(genes) if genes else features_count.get("gene", 0)

    return {
        "n_genes": n_genes,
        "n_transcripts": n_transcripts,
        "n_exons": features_count.get("exon", 0),
        "n_cds": features_count.get("CDS", 0),
        "mono_exonic_genes": mono_exonic,
        "mean_exons_per_transcript": mean(exon_counts) if exon_counts else 0,
        "median_exons_per_transcript": median(exon_counts) if exon_counts else 0,
        "intron_lengths": intron_lengths,
    }


def parse_busco(path: Path) -> Dict[str, float]:
    values = {
        "busco_complete": 0.0,
        "busco_single_copy": 0.0,
        "busco_duplicated": 0.0,
        "busco_fragmented": 0.0,
        "busco_missing": 0.0,
        "busco_n_markers": 0.0,
    }
    try:
        data = json.loads(path.read_text())
    except Exception:
        return values

    # BUSCO JSON structures vary between versions; support common layouts.
    results = data.get("results", data)
    one_line = results.get("one_line_summary_raw", data.get("one_line_summary_raw", ""))
    if isinstance(one_line, str):
        for key, outkey in [("C", "busco_complete"), ("S", "busco_single_copy"), ("D", "busco_duplicated"), ("F", "busco_fragmented"), ("M", "busco_missing"), ("n", "busco_n_markers")]:
            m = re.search(rf"{key}:([0-9.]+)%?", one_line)
            if m:
                values[outkey] = float(m.group(1))

    for src in (data, results):
        for k, outkey in [
            ("Complete percentage", "busco_complete"),
            ("Single copy percentage", "busco_single_copy"),
            ("Multi copy percentage", "busco_duplicated"),
            ("Fragmented percentage", "busco_fragmented"),
            ("Missing percentage", "busco_missing"),
            ("n_markers", "busco_n_markers"),
        ]:
            if k in src:
                try:
                    values[outkey] = float(src[k])
                except Exception:
                    pass
    return values


def parse_omark(path: Path) -> Dict[str, float]:
    values = {
        "omark_conserved": 0.0,
        "omark_single": 0.0,
        "omark_duplicated": 0.0,
        "omark_fragmented": 0.0,
        "omark_missing": 0.0,
        "omark_contaminant": 0.0,
        "omark_unknown": 0.0,
    }
    text = path.read_text(errors="ignore") if path.exists() else ""
    patterns = {
        "omark_conserved": r"conserved[^0-9]*([0-9.]+)",
        "omark_single": r"single[^0-9]*([0-9.]+)",
        "omark_duplicated": r"duplicated[^0-9]*([0-9.]+)",
        "omark_fragmented": r"fragmented[^0-9]*([0-9.]+)",
        "omark_missing": r"missing[^0-9]*([0-9.]+)",
        "omark_contaminant": r"contamin[^0-9]*([0-9.]+)",
        "omark_unknown": r"unknown[^0-9]*([0-9.]+)",
    }
    for key, pat in patterns.items():
        m = re.search(pat, text, flags=re.IGNORECASE)
        if m:
            try:
                values[key] = float(m.group(1))
            except Exception:
                pass
    return values


def write_summary(prefix: str, protein_lengths: Dict[str, int], gff_stats: dict, busco: dict, omark: dict) -> None:
    lengths = list(protein_lengths.values())
    n_prot = len(lengths)
    short = sum(1 for x in lengths if x < 100)
    row = {
        "genome_id": prefix.split("_")[0],
        "annotation_method": "_".join(prefix.split("_")[1:-1]) if "_" in prefix else "NA",
        "n_genes": gff_stats["n_genes"],
        "n_transcripts": gff_stats["n_transcripts"],
        "n_proteins": n_prot,
        "mean_protein_length": round(mean(lengths), 3) if lengths else 0,
        "median_protein_length": round(median(lengths), 3) if lengths else 0,
        "short_proteins_lt100": short,
        "short_proteins_lt100_percent": round(100 * short / n_prot, 3) if n_prot else 0,
        "mono_exonic_genes": gff_stats["mono_exonic_genes"],
        "mean_exons_per_transcript": round(gff_stats["mean_exons_per_transcript"], 3),
        "median_exons_per_transcript": round(gff_stats["median_exons_per_transcript"], 3),
        **busco,
        **omark,
    }
    out = Path(f"{prefix}_summary.tsv")
    with out.open("w") as handle:
        handle.write("\t".join(row.keys()) + "\n")
        handle.write("\t".join(map(str, row.values())) + "\n")


def write_distributions(prefix: str, protein_lengths: Dict[str, int], intron_lengths: List[int]) -> None:
    with Path(f"{prefix}_protein_lengths.tsv").open("w") as handle:
        handle.write("protein_id\tprotein_length\n")
        for pid, length in sorted(protein_lengths.items()):
            handle.write(f"{pid}\t{length}\n")
    with Path(f"{prefix}_intron_lengths.tsv").open("w") as handle:
        handle.write("intron_id\tintron_length\n")
        for i, length in enumerate(intron_lengths, start=1):
            handle.write(f"intron_{i}\t{length}\n")


def main() -> None:
    args = parse_args()
    proteins = fasta_lengths(Path(args.proteins))
    gff_stats = parse_gff(Path(args.gff))
    busco = parse_busco(Path(args.busco))
    omark = parse_omark(Path(args.omark))
    write_summary(args.prefix, proteins, gff_stats, busco, omark)
    write_distributions(args.prefix, proteins, gff_stats["intron_lengths"])


if __name__ == "__main__":
    main()
