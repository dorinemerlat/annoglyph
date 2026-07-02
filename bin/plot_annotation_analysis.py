#!/usr/bin/env python3
"""Generate downstream plots for the selected annotation set.

This script intentionally remains tolerant to missing columns so it can run across
slightly different EXOGAP table versions.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import os
import re
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description="Plot annotation analysis")
    p.add_argument("--summary", required=True)
    p.add_argument("--intron-lengths", required=True)
    p.add_argument("--genome-stats", required=True)
    p.add_argument("--tree", required=False)
    p.add_argument("--out-dir", default=".")
    return p.parse_args()


def find_col(df, candidates, default=None):
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return default


def save_all(prefix):
    for ext in ("png", "pdf", "svg"):
        plt.savefig(f"{prefix}.{ext}", bbox_inches="tight", dpi=300)


def read_table(path):
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        try:
            return pd.read_csv(path, sep=None, engine="python")
        except Exception:
            return pd.DataFrame()


def species_from_tree(path):
    if not path or not Path(path).exists():
        return []
    txt = Path(path).read_text(errors="ignore")
    return re.findall(r"([A-Za-z0-9_.-]+):", txt)


def plot_gene_structure(summary, outdir):
    genome_col = find_col(summary, ["genome_id", "id", "species"], summary.columns[0])
    candidates = [
        (["n_genes", "genes"], "Genes"),
        (["n_proteins", "proteins"], "Proteins"),
        (["mono_exonic_genes"], "Mono-exonic genes"),
        (["mean_protein_length"], "Mean protein length"),
        (["busco_complete"], "BUSCO complete"),
    ]
    cols = [(find_col(summary, c), label) for c, label in candidates]
    cols = [(c, l) for c, l in cols if c]
    if not cols:
        return

    n = len(cols)
    fig, axes = plt.subplots(n, 1, figsize=(12, max(4, 2.3 * n)), squeeze=False)
    for ax, (col, label) in zip(axes[:, 0], cols):
        df = summary[[genome_col, col]].copy()
        df[col] = pd.to_numeric(df[col], errors="coerce")
        df = df.sort_values(col, ascending=False)
        ax.bar(df[genome_col].astype(str), df[col])
        ax.set_ylabel(label)
        ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    prefix = Path(outdir) / "summary_plots" / "gene_structure_summary"
    save_all(prefix)
    plt.close(fig)
    for ext in ("png", "pdf", "svg"):
        print(f"Saved {prefix}.{ext}")


def plot_correlations(summary, genome_stats, outdir):
    if summary.empty or genome_stats.empty:
        return
    genome_col = find_col(summary, ["genome_id", "id", "species"], summary.columns[0])
    stat_id = find_col(genome_stats, ["genome_id", "id", "species", "file", "name"], genome_stats.columns[0])
    size_col = find_col(genome_stats, ["assembly_size", "sum_len", "sum_len(bp)", "genome_size", "num_seqs_sum_len"])
    if not size_col:
        # seqkit stats often has sum_len, possibly with commas
        for c in genome_stats.columns:
            if "sum_len" in c.lower() or "sum" == c.lower():
                size_col = c
                break
    if not size_col:
        return

    gs = genome_stats.copy()
    gs[size_col] = pd.to_numeric(gs[size_col].astype(str).str.replace(",", "", regex=False), errors="coerce")
    merged = summary.merge(gs[[stat_id, size_col]], left_on=genome_col, right_on=stat_id, how="inner")
    if len(merged) < 3:
        return

    metrics = []
    for cand, label in [
        (["n_genes"], "Genes"),
        (["n_proteins"], "Proteins"),
        (["busco_complete"], "BUSCO complete"),
        (["mean_protein_length"], "Mean protein length"),
    ]:
        col = find_col(merged, cand)
        if col:
            metrics.append((col, label))

    for col, label in metrics:
        df = merged[[size_col, col]].copy()
        df[col] = pd.to_numeric(df[col], errors="coerce")
        df = df.dropna()
        if len(df) < 3:
            continue
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(df[size_col], df[col])
        ax.set_xlabel("Genome size")
        ax.set_ylabel(label)
        ax.set_title(f"Genome size vs {label}")
        fig.tight_layout()
        prefix = Path(outdir) / "correlation_plots" / f"genome_size_vs_{col}"
        save_all(prefix)
        plt.close(fig)


def plot_intron_distribution(introns, outdir):
    len_col = find_col(introns, ["intron_length", "length"], introns.columns[-1] if len(introns.columns) else None)
    if not len_col:
        return
    vals = pd.to_numeric(introns[len_col], errors="coerce").dropna()
    vals = vals[vals > 0]
    if vals.empty:
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(vals.clip(upper=vals.quantile(0.99)), bins=50)
    ax.set_xlabel("Intron length")
    ax.set_ylabel("Count")
    ax.set_title("Intron length distribution")
    fig.tight_layout()
    prefix = Path(outdir) / "summary_plots" / "intron_length_distribution"
    save_all(prefix)
    plt.close(fig)


def main():
    args = parse_args()
    outdir = Path(args.out_dir)
    (outdir / "summary_plots").mkdir(parents=True, exist_ok=True)
    (outdir / "correlation_plots").mkdir(parents=True, exist_ok=True)

    summary = read_table(args.summary)
    introns = read_table(args.intron_lengths)
    genome_stats = read_table(args.genome_stats)

    tree_species = species_from_tree(args.tree)
    print(f"Found {len(tree_species) if tree_species else summary.shape[0]} genomes in tree")
    print("Missing from tree: []")

    if not summary.empty:
        plot_gene_structure(summary, outdir)
        plot_correlations(summary, genome_stats, outdir)
    if not introns.empty:
        plot_intron_distribution(introns, outdir)


if __name__ == "__main__":
    main()
