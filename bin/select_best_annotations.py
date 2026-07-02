#!/usr/bin/env python3
"""Select the best annotation method per genome from a merged summary table."""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Select best annotation per genome")
    p.add_argument("--summary", required=True, help="Merged summary TSV")
    p.add_argument("--output", "-o", default="best_annotations.tsv", help="Output TSV")
    return p.parse_args()


def find_col(df, candidates, default=None):
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return default


def to_num(s):
    return pd.to_numeric(s, errors="coerce").fillna(0)


def main():
    args = parse_args()
    df = pd.read_csv(args.summary, sep="\t")
    if df.empty:
        raise SystemExit("ERROR: summary table is empty")

    genome_col = find_col(df, ["genome_id", "id", "species", "species_id"], df.columns[0])
    method_col = find_col(df, ["annotation_method", "method", "annotation"], df.columns[1] if len(df.columns) > 1 else df.columns[0])

    busco_col = find_col(df, ["busco_complete", "complete_busco", "Complete BUSCOs (%)"])
    missing_col = find_col(df, ["busco_missing", "missing_busco", "Missing BUSCOs (%)"])
    omark_cont_col = find_col(df, ["omark_contaminant", "omark_contamination", "contaminant"])
    nprot_col = find_col(df, ["n_proteins", "proteins", "protein_count"])
    short_col = find_col(df, ["short_proteins_lt100_percent", "short_proteins_percent"])
    mono_col = find_col(df, ["mono_exonic_genes", "monoexonic_genes"])

    score = pd.Series(0.0, index=df.index)
    if busco_col:
        score += to_num(df[busco_col]) * 2.0
    if missing_col:
        score -= to_num(df[missing_col]) * 1.0
    if omark_cont_col:
        score -= to_num(df[omark_cont_col]) * 2.0
    if short_col:
        score -= to_num(df[short_col]) * 0.5
    if mono_col:
        # Penalize extreme mono-exonic inflation gently after scaling.
        mono = to_num(df[mono_col])
        if mono.max() > 0:
            score -= 5 * mono / mono.max()
    if nprot_col:
        nprot = to_num(df[nprot_col])
        median = nprot.median() if len(nprot) else 0
        if median > 0:
            score -= ((nprot - median).abs() / median) * 5

    df = df.copy()
    df["selection_score"] = score

    rows = []
    for genome_id, sub in df.groupby(genome_col, sort=False):
        best = sub.sort_values("selection_score", ascending=False).iloc[0]
        method = best[method_col]
        sufficient_orthoinspector = "yes"
        sufficient_orthofinder = "yes"
        if busco_col and float(best.get(busco_col, 0)) < 50:
            sufficient_orthoinspector = "no"
            sufficient_orthofinder = "no"
        rows.append({
            "genome_id": genome_id,
            "annotation_method": method,
            "sufficient_quality_for_orthoinspector": sufficient_orthoinspector,
            "sufficient_quality_for_orthofinder": sufficient_orthofinder,
            "selection_score": round(float(best["selection_score"]), 6),
        })

    out = pd.DataFrame(rows)
    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
