#!/usr/bin/env python3
"""Plot compact summaries comparing gene annotation methods.

Expected inputs:
  --summary merged summary TSV
  --protein-lengths merged protein length TSV
  --tree optional Newick tree
  --out output prefix
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description="Plot annotation method comparison")
    p.add_argument("--summary", required=True)
    p.add_argument("--protein-lengths", required=True)
    p.add_argument("--tree", required=False)
    p.add_argument("--out", required=True)
    return p.parse_args()


def find_col(df, candidates, default=None):
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return default


def species_from_tree(path):
    if not path or not Path(path).exists():
        return []
    txt = Path(path).read_text(errors="ignore")
    names = re.findall(r"([A-Za-z0-9_.-]+):", txt)
    return names


def save_all(prefix):
    for ext in ("png", "pdf", "svg"):
        plt.savefig(f"{prefix}.{ext}", bbox_inches="tight", dpi=300)


def main():
    args = parse_args()
    summary = pd.read_csv(args.summary, sep="\t")
    prot = pd.read_csv(args.protein_lengths, sep="\t")

    if summary.empty:
        raise SystemExit("ERROR: empty summary table")

    genome_col = find_col(summary, ["genome_id", "id", "species"], summary.columns[0])
    method_col = find_col(summary, ["annotation_method", "method"], summary.columns[1] if len(summary.columns) > 1 else summary.columns[0])

    tree_species = species_from_tree(args.tree)
    missing = []
    if tree_species:
        available = set(map(str, summary[genome_col]))
        missing = [x for x in tree_species if x not in available]
    print(f"Found {len(tree_species) if tree_species else summary[genome_col].nunique()} genomes in tree")
    print(f"Missing from tree: {missing}")

    metrics = []
    for cands, label in [
        (["n_proteins", "proteins", "protein_count"], "Proteins"),
        (["n_genes", "genes", "gene_count"], "Genes"),
        (["busco_complete", "Complete BUSCOs (%)"], "BUSCO complete"),
        (["omark_contaminant", "omark_contamination"], "OMArk contaminants"),
        (["short_proteins_lt100_percent"], "Short proteins (%)"),
    ]:
        col = find_col(summary, cands)
        if col:
            tmp = summary[[genome_col, method_col, col]].copy()
            tmp["metric"] = label
            tmp = tmp.rename(columns={col: "value"})
            tmp["value"] = pd.to_numeric(tmp["value"], errors="coerce")
            metrics.append(tmp)

    if not metrics:
        raise SystemExit("ERROR: no recognized metrics in summary table")

    mdf = pd.concat(metrics, ignore_index=True)
    agg = mdf.groupby(["metric", method_col], dropna=False)["value"].median().reset_index()

    n_metrics = agg["metric"].nunique()
    fig_h = max(4, 2.4 * n_metrics)
    fig, axes = plt.subplots(n_metrics, 1, figsize=(10, fig_h), squeeze=False)

    for ax, (metric, sub) in zip(axes[:, 0], agg.groupby("metric", sort=False)):
        sub = sub.sort_values("value", ascending=False)
        ax.bar(sub[method_col].astype(str), sub["value"])
        ax.set_title(metric)
        ax.set_ylabel("Median")
        ax.tick_params(axis="x", rotation=35)

    fig.suptitle("Annotation comparison summary", y=1.01)
    fig.tight_layout()
    save_all(args.out)
    plt.close(fig)
    print(f"Saved {args.out}.png")
    print(f"Saved {args.out}.pdf")
    print(f"Saved {args.out}.svg")


if __name__ == "__main__":
    main()
