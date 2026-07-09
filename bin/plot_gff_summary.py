#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns

try:
    from Bio import Phylo
except ImportError:
    Phylo = None


def clean_genome_name(genome_id: str) -> str:
    name = re.sub(r"-\d+$", "", str(genome_id))
    name = name.replace("-", " ")
    parts = name.split()
    if parts:
        parts[0] = parts[0].capitalize()
    return " ".join(parts)


def read_table(tables_dir: Path, name: str) -> pd.DataFrame:
    path = tables_dir / name
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def get_genome_order(all_genomes: list[str], tree_file: Path | None) -> list[str]:
    if tree_file is None:
        return sorted(all_genomes)

    if Phylo is None:
        raise ImportError("Biopython required: pip install biopython")

    tree = Phylo.read(tree_file, "newick")
    phylo_order = [leaf.name for leaf in tree.get_terminals()]
    phylo_order = [g for g in phylo_order if g in all_genomes]
    missing = sorted([g for g in all_genomes if g not in phylo_order])

    return phylo_order + missing


def setup_xlabels(ax, genome_order: list[str], labelsize: int = 6) -> None:
    ax.set_xticks(range(len(genome_order)))
    ax.set_xticklabels(
        [clean_genome_name(g) for g in genome_order],
        rotation=90,
        fontsize=labelsize,
        fontstyle="italic",
    )


def apply_order(df: pd.DataFrame, genome_order: list[str]) -> pd.DataFrame:
    if df.empty:
        return df

    df = df.copy()
    df["genome_id"] = pd.Categorical(
        df["genome_id"],
        categories=genome_order,
        ordered=True,
    )
    return df


def make_gff_figure(
    method: str,
    genome_order: list[str],
    out_dir: Path,
    df_gff: pd.DataFrame,
) -> None:

    if df_gff.empty:
        print("No gff_stats.tsv found, skipping.")
        return

    start_cols = [
        "percent_genes_start_and_stop",
        "percent_genes_start_only",
        "percent_genes_stop_only",
    ]

    annot_cols = [
        "percent_genes_product_only",
        "percent_genes_go_only",
        "percent_genes_product_and_go",
    ]

    required = (
        start_cols
        + [
            "mean_exons_per_gene",
            "percent_monoexonic_genes",
        ]
        + annot_cols
    )

    missing = [c for c in required if c not in df_gff.columns]

    if missing:
        print(f"Skipping GFF figure for {method}: missing columns {missing}")
        return

    df_m = (
        df_gff[df_gff["method_annotation"] == method]
        .groupby("genome_id", as_index=False, observed=False)[required]
        .mean()
        .set_index("genome_id")
        .reindex(genome_order)
        .reset_index()
    )

    df_m = apply_order(df_m, genome_order)

    fig, axes = plt.subplots(
        nrows=4,
        ncols=1,
        figsize=(8.27, 11.69),
        sharex=True,
        gridspec_kw={"height_ratios": [1.2, 1, 1, 1.2]},
    )

    df_m[["genome_id"] + start_cols].set_index("genome_id").plot(
        kind="bar",
        stacked=True,
        ax=axes[0],
        width=0.9,
    )
    axes[0].set_title("Start and stop codon annotation")
    axes[0].set_ylabel("Percentage (%)")
    axes[0].set_ylim(0, 100)
    axes[0].yaxis.set_major_formatter(mticker.PercentFormatter())
    axes[0].legend(
        ["Start + stop", "Start only", "Stop only"],
        title="Category",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=7,
    )

    sns.barplot(
        data=df_m,
        x="genome_id",
        y="mean_exons_per_gene",
        order=genome_order,
        ax=axes[1],
        color="#8c8c8c",
    )
    axes[1].set_title("Mean number of exons per gene")
    axes[1].set_ylabel("Mean exons/gene")

    sns.barplot(
        data=df_m,
        x="genome_id",
        y="percent_monoexonic_genes",
        order=genome_order,
        ax=axes[2],
        color="#8c8c8c",
    )
    axes[2].set_title("Mono-exonic genes")
    axes[2].set_ylabel("Percentage (%)")
    axes[2].set_ylim(0, 100)
    axes[2].yaxis.set_major_formatter(mticker.PercentFormatter())

    df_m[["genome_id"] + annot_cols].set_index("genome_id").plot(
        kind="bar",
        stacked=True,
        ax=axes[3],
        width=0.9,
    )
    axes[3].set_title("Functional annotation")
    axes[3].set_ylabel("Percentage (%)")
    axes[3].set_ylim(0, 100)
    axes[3].yaxis.set_major_formatter(mticker.PercentFormatter())
    axes[3].legend(
        ["Product only", "GO only", "Product + GO"],
        title="Category",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=7,
    )

    for ax in axes:
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelsize=6)

    setup_xlabels(axes[-1], genome_order)
    axes[-1].set_xlabel("Genome ID")

    if method == "all_annotations":
        fig.suptitle("GFF metrics", y=1.005, fontsize=12)
    else:
        fig.suptitle(f"{method} — GFF metrics", y=1.005, fontsize=12)

    plt.tight_layout()

    out = out_dir / f"{method}_gff_summary.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {out}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tables-dir", required=True, type=Path)
    parser.add_argument("-o", "--out-dir", required=True, type=Path)
    parser.add_argument("--tree", type=Path, default=None)
    parser.add_argument("--methods", nargs="*", default=None)
    parser.add_argument("--ignore-method", action="store_true")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    df_gff = read_table(args.tables_dir, "gff_stats.tsv")

    if df_gff.empty:
        raise FileNotFoundError(f"No gff_stats.tsv found in {args.tables_dir}")

    if args.ignore_method and "method_annotation" in df_gff.columns:
        df_gff["method_annotation"] = "all_annotations"

    all_genomes = sorted(set(df_gff.get("genome_id", [])))
    genome_order = get_genome_order(all_genomes, args.tree)

    if args.ignore_method:
        methods = ["all_annotations"]
    elif args.methods:
        methods = args.methods
    else:
        methods = sorted(set(df_gff.get("method_annotation", [])))

    for method in methods:
        make_gff_figure(
            method,
            genome_order,
            args.out_dir,
            df_gff,
        )


if __name__ == "__main__":
    main()