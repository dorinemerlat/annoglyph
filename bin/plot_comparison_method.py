#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.ticker as mticker

try:
    from Bio import Phylo
except ImportError:
    Phylo = None


METHOD_ORDER = [
    "galba",
    "braker2",
    "braker_with_varus",
    "braker_with_star",
]

METHOD_COLORS = {
    "galba": "#1b9e77",
    "braker2": "#d95f02",
    "braker_with_varus": "#e7298a",
    "braker_with_star": "#7570b3",
}

METHOD_MARKERS = {
    "galba": "o",
    "braker2": "^",
    "braker_with_varus": "D",
    "braker_with_star": "s",
}

NEUTRAL_BAR_COLOR = "#8c8c8c"
NEUTRAL_BOX_COLOR = "#bdbdbd"

BUSCO_COLORS = {
    "single": "#69c4ec",
    "duplicated": "#3fa7d6",
    "fragmented": "#f0e442",
    "missing": "#f45b55",
}

LABEL_SIZE = 8
TITLE_SIZE = 10
MAIN_TITLE_SIZE = 13
LEGEND_SIZE = 8
LEGEND_TITLE_SIZE = 8

plt.rcParams["hatch.linewidth"] = 0.4


def clean_genome_name(genome_id: str) -> str:
    name = re.sub(r"-\d$", "", str(genome_id))
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
        raise ImportError(
            "Biopython is required to read Newick files. Install with: pip install biopython"
        )

    tree = Phylo.read(tree_file, "newick")

    phylo_order = []
    for leaf in tree.get_terminals():
        genome_id = str(leaf.name).split("|")[0]
        if genome_id in all_genomes:
            phylo_order.append(genome_id)

    missing = sorted([g for g in all_genomes if g not in phylo_order])

    print(f"Found {len(phylo_order)} genomes in tree")
    print(f"Missing from tree: {missing}")

    return phylo_order + missing


def get_methods(*dfs: pd.DataFrame) -> list[str]:
    methods = set()
    for df in dfs:
        if not df.empty and "method_annotation" in df.columns:
            methods.update(df["method_annotation"].dropna().astype(str))
    return sorted(methods)


def build_positions(
    genome_order: list[str],
    methods: list[str],
    available_pairs: set[tuple[str, str]],
) -> tuple[pd.DataFrame, dict[str, float]]:
    rows = []
    genome_centers = {}

    x = 0.0
    bar_spacing = 0.22
    genome_gap = 0.35

    for genome_id in genome_order:
        available_methods = [
            method for method in methods
            if (genome_id, method) in available_pairs
        ]

        if not available_methods:
            continue

        xs = []

        for method in available_methods:
            rows.append({
                "genome_id": genome_id,
                "method_annotation": method,
                "x": x,
            })
            xs.append(x)
            x += bar_spacing

        genome_centers[genome_id] = sum(xs) / len(xs)
        x += genome_gap

    return pd.DataFrame(rows), genome_centers


def add_method_markers(ax, pos_df: pd.DataFrame, y: float = -8) -> None:
    for _, row in pos_df.iterrows():
        method = row["method_annotation"]
        ax.scatter(
            row["x"],
            y,
            marker=METHOD_MARKERS.get(method, "o"),
            color=METHOD_COLORS.get(method, "gray"),
            s=18,
            clip_on=False,
            zorder=10,
        )


def setup_genome_xaxis(ax, genome_centers: dict[str, float]) -> None:
    genomes = list(genome_centers.keys())

    ax.set_xticks([genome_centers[g] for g in genomes])
    ax.set_xticklabels(
        [clean_genome_name(g) for g in genomes],
        rotation=90,
        fontsize=LABEL_SIZE,
        fontstyle="italic",
    )


def method_legend_handles(methods: list[str]) -> list[mlines.Line2D]:
    return [
        mlines.Line2D(
            [],
            [],
            color=METHOD_COLORS.get(method, "gray"),
            marker=METHOD_MARKERS.get(method, "o"),
            linestyle="None",
            markersize=6,
            label=method,
        )
        for method in methods
    ]


def style_axis(ax) -> None:
    ax.tick_params(axis="both", labelsize=LABEL_SIZE)
    ax.xaxis.label.set_size(LABEL_SIZE)
    ax.yaxis.label.set_size(LABEL_SIZE)
    ax.xaxis.label.set_weight("bold")
    ax.yaxis.label.set_weight("bold")
    ax.title.set_size(TITLE_SIZE)


def style_legend(legend) -> None:
    if legend is None:
        return

    legend.get_title().set_fontsize(LEGEND_TITLE_SIZE)
    legend.get_title().set_fontweight("bold")

    for text in legend.get_texts():
        text.set_fontsize(LEGEND_SIZE)


def plot_grouped_bar(
    ax,
    df: pd.DataFrame,
    value_col: str,
    pos_df: pd.DataFrame,
    title: str,
    ylabel: str,
    marker_y: float | None = None,
) -> None:
    plot_df = pos_df.merge(df, on=["genome_id", "method_annotation"], how="left")

    ax.bar(
        plot_df["x"],
        plot_df[value_col],
        width=0.18,
        color=NEUTRAL_BAR_COLOR,
        edgecolor="none",
    )

    ax.set_title(title, fontsize=TITLE_SIZE, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=LABEL_SIZE, fontweight="bold")
    ax.set_xlabel("")

    if marker_y is not None:
        add_method_markers(ax, plot_df, y=marker_y)


def plot_grouped_boxplot(
    ax,
    df: pd.DataFrame,
    pos_df: pd.DataFrame,
) -> None:
    data = []
    positions = []

    for _, row in pos_df.iterrows():
        genome_id = row["genome_id"]
        method = row["method_annotation"]

        values = df[
            (df["genome_id"] == genome_id)
            & (df["method_annotation"] == method)
        ]["protein_length"].dropna()

        if len(values) == 0:
            continue

        data.append(values)
        positions.append(row["x"])

    bp = ax.boxplot(
        data,
        positions=positions,
        widths=0.16,
        showfliers=False,
        patch_artist=True,
    )

    for patch in bp["boxes"]:
        patch.set_facecolor(NEUTRAL_BOX_COLOR)
        patch.set_alpha(0.8)

    for element in ["whiskers", "caps", "medians"]:
        for item in bp[element]:
            item.set_color("black")
            item.set_linewidth(0.6)

    ax.set_title("Protein length distribution", fontsize=TITLE_SIZE, fontweight="bold")
    ax.set_ylabel("Protein length (aa)", fontsize=LABEL_SIZE, fontweight="bold")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f"{int(x):,}"))

    add_method_markers(ax, pos_df, y=-120)


def plot_stacked_grouped_bar(
    ax,
    df: pd.DataFrame,
    pos_df: pd.DataFrame,
    value_cols: list[str],
    colors: list[str],
    labels: list[str],
    title: str,
    ylabel: str = "Percentage (%)",
    hatches: list[str] | None = None,
    edges: list[str] | None = None,
) -> None:
    plot_df = pos_df.merge(df, on=["genome_id", "method_annotation"], how="left")

    bottoms = pd.Series([0.0] * len(plot_df))

    for i, col in enumerate(value_cols):
        values = plot_df[col].fillna(0)

        bars = ax.bar(
            plot_df["x"],
            values,
            width=0.18,
            bottom=bottoms,
            color=colors[i],
            label=labels[i],
        )

        if hatches is not None:
            for bar in bars:
                hatch = hatches[i]
                if hatch:
                    bar.set_hatch(hatch)
                    bar.set_edgecolor(edges[i] if edges is not None else "black")
                    bar.set_linewidth(0)
                else:
                    bar.set_edgecolor("none")
                    bar.set_linewidth(0)

        bottoms += values

    ax.set_title(title, fontsize=TITLE_SIZE, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=LABEL_SIZE, fontweight="bold")
    ax.set_ylim(-12, 100)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter())

    add_method_markers(ax, plot_df, y=-7)


def make_combined_figure(
    genome_order: list[str],
    methods: list[str],
    available_pairs: set[tuple[str, str]],
    out_dir: Path,
    df_gene_number: pd.DataFrame,
    df_protein_lengths: pd.DataFrame,
    df_busco: pd.DataFrame,
    df_omark: pd.DataFrame,
    df_gff: pd.DataFrame,
) -> None:
    pos_df, genome_centers = build_positions(genome_order, methods, available_pairs)

    busco_cols = [
        "busco_single_percent",
        "busco_duplicated_percent",
        "busco_fragmented_percent",
        "busco_missing_percent",
    ]

    busco_labels = [
        "Single-copy",
        "Duplicated",
        "Fragmented",
        "Missing",
    ]

    busco_colors = [
        BUSCO_COLORS["single"],
        BUSCO_COLORS["duplicated"],
        BUSCO_COLORS["fragmented"],
        BUSCO_COLORS["missing"],
    ]

    omark_hog_cols = [
        "omark_single_percent",
        "omark_duplicated_percent",
        "omark_missing_percent",
    ]

    omark_hog_labels = [
        "Single",
        "Duplicated",
        "Missing",
    ]

    omark_hog_colors = [
        "#46bea1ff",
        "#cfee8eff",
        "#ed1c5aff",
    ]

    omark_consistency_cols = [
        "omark_consistent_complete_percent",
        "omark_consistent_partial_percent",
        "omark_consistent_fragmented_percent",
        "omark_contamination_complete_percent",
        "omark_contamination_partial_percent",
        "omark_contamination_fragmented_percent",
        "omark_inconsistent_complete_percent",
        "omark_inconsistent_partial_percent",
        "omark_inconsistent_fragmented_percent",
        "omark_unknown_percent",
    ]

    omark_consistency_colors = [
        "#3db7e9ff",
        "#3db7e9ff",
        "#3db7e9ff",
        "#e69f00ff",
        "#e69f00ff",
        "#e69f00ff",
        "#8a5df1ff",
        "#8a5df1ff",
        "#8a5df1ff",
        "#000000ff",
    ]

    omark_hatches = [
        "",
        "///",
        "///",
        "",
        "///",
        "///",
        "",
        "///",
        "///",
        "",
    ]

    omark_edges = [
        "none",
        "darkslategray",
        "white",
        "none",
        "darkslategray",
        "white",
        "none",
        "darkslategray",
        "white",
        "none",
    ]

    df_genes = (
        df_gene_number
        .groupby(["genome_id", "method_annotation"], as_index=False, observed=False)["gene_number"]
        .sum()
    )

    df_lengths = df_protein_lengths.copy()

    df_short = (
        df_lengths
        .assign(is_short=lambda d: d["protein_length"] < 100)
        .groupby(["genome_id", "method_annotation"], observed=False)["is_short"]
        .mean()
        .mul(100)
        .reset_index(name="percent_proteins_lt_100aa")
    )

    if (
        not df_gff.empty
        and {"percent_genes_start_and_stop", "percent_genes_start_only"}.issubset(df_gff.columns)
    ):
        df_start = (
            df_gff
            .assign(
                percent_genes_with_start=lambda d: (
                    d["percent_genes_start_and_stop"].fillna(0)
                    + d["percent_genes_start_only"].fillna(0)
                )
            )
            .groupby(["genome_id", "method_annotation"], as_index=False, observed=False)["percent_genes_with_start"]
            .mean()
        )
    else:
        df_start = pd.DataFrame(
            columns=["genome_id", "method_annotation", "percent_genes_with_start"]
        )

    df_busco_m = (
        df_busco
        .groupby(["genome_id", "method_annotation"], as_index=False, observed=False)[busco_cols]
        .mean()
    )

    df_omark_m = (
        df_omark
        .groupby(["genome_id", "method_annotation"], as_index=False, observed=False)[
            omark_hog_cols + omark_consistency_cols
        ]
        .mean()
    )

    fig, axes = plt.subplots(
        nrows=7,
        ncols=1,
        figsize=(22, 11.69),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 1, 1, 1, 1, 1, 1]},
    )

    plot_grouped_bar(
        axes[0],
        df_genes,
        "gene_number",
        pos_df,
        "Number of predicted genes",
        "Gene count",
        marker_y=-7000,
    )
    axes[0].yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f"{int(x):,}"))

    plot_grouped_boxplot(
        axes[1],
        df_lengths,
        pos_df,
    )

    plot_grouped_bar(
        axes[2],
        df_short,
        "percent_proteins_lt_100aa",
        pos_df,
        "Short proteins (<100 aa)",
        "Percentage (%)",
        marker_y=-7,
    )
    axes[2].yaxis.set_major_formatter(mticker.PercentFormatter())
    axes[2].set_ylim(-12, max(axes[2].get_ylim()[1], 60))

    plot_grouped_bar(
        axes[3],
        df_start,
        "percent_genes_with_start",
        pos_df,
        "Genes with a codon start annotation",
        "Percentage (%)",
        marker_y=-7,
    )
    axes[3].yaxis.set_major_formatter(mticker.PercentFormatter())
    axes[3].set_ylim(-12, 100)

    plot_stacked_grouped_bar(
        axes[4],
        df_busco_m,
        pos_df,
        busco_cols,
        busco_colors,
        busco_labels,
        "BUSCO completeness",
    )

    busco_legend = axes[4].legend(
        handles=[
            mpatches.Patch(facecolor=color, label=label)
            for color, label in zip(busco_colors, busco_labels)
        ],
        title="BUSCO category",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=LEGEND_SIZE,
        title_fontsize=LEGEND_TITLE_SIZE,
    )
    style_legend(busco_legend)

    plot_stacked_grouped_bar(
        axes[5],
        df_omark_m,
        pos_df,
        omark_hog_cols,
        omark_hog_colors,
        omark_hog_labels,
        "OMArk conserved HOGs",
    )

    omark_hog_legend = axes[5].legend(
        handles=[
            mpatches.Patch(facecolor=color, label=label)
            for color, label in zip(omark_hog_colors, omark_hog_labels)
        ],
        title="Completeness",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=LEGEND_SIZE,
        title_fontsize=LEGEND_TITLE_SIZE,
    )
    style_legend(omark_hog_legend)

    plot_stacked_grouped_bar(
        axes[6],
        df_omark_m,
        pos_df,
        omark_consistency_cols,
        omark_consistency_colors,
        [
            "Consistent complete",
            "Consistent partial",
            "Consistent fragmented",
            "Contamination complete",
            "Contamination partial",
            "Contamination fragmented",
            "Inconsistent complete",
            "Inconsistent partial",
            "Inconsistent fragmented",
            "Unknown",
        ],
        "OMArk proteome consistency",
        hatches=omark_hatches,
        edges=omark_edges,
    )

    consistency_legend = axes[6].legend(
        handles=[
            mpatches.Patch(facecolor="#3db7e9ff", edgecolor="white", label="Consistent"),
            mpatches.Patch(facecolor="#e69f00ff", edgecolor="white", label="Contamination"),
            mpatches.Patch(facecolor="#8a5df1ff", edgecolor="white", label="Inconsistent"),
            mpatches.Patch(facecolor="#000000ff", edgecolor="white", label="Unknown"),
            mpatches.Patch(facecolor="white", hatch="///", edgecolor="darkslategray", label="Partial mapping"),
            mpatches.Patch(facecolor="darkslategray", hatch="///", edgecolor="white", label="Fragments"),
        ],
        title="Consistency",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=LEGEND_SIZE,
        title_fontsize=LEGEND_TITLE_SIZE,
    )
    style_legend(consistency_legend)

    x_min = pos_df["x"].min() - 0.3
    x_max = pos_df["x"].max() + 0.3

    for ax in axes:
        ax.set_xlim(x_min, x_max)
        ax.set_xlabel("")
        style_axis(ax)

    # Align all Y-axis labels vertically
    for ax in axes:
        ax.yaxis.set_label_coords(-0.040, 0.5)

    setup_genome_xaxis(axes[-1], genome_centers)
    axes[-1].set_xlabel("Genome ID", fontsize=LABEL_SIZE, fontweight="bold")

    method_legend = fig.legend(
        handles=method_legend_handles(methods),
        title="Annotation method",
        loc="lower center",
        bbox_to_anchor=(0.5, 0.01),
        ncol=min(len(methods), 4),
        fontsize=LEGEND_SIZE,
        title_fontsize=LEGEND_TITLE_SIZE,
    )
    style_legend(method_legend)

    plt.tight_layout(rect=[0.05, 0.06, 0.98, 0.99])

    out = out_dir / "combined_protein_busco_omark_summary.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {out}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate one combined protein/BUSCO/OMArk summary figure."
    )

    parser.add_argument(
        "-t",
        "--tables-dir",
        required=True,
        type=Path,
        help="Directory containing TSV tables.",
    )

    parser.add_argument(
        "-o",
        "--out-dir",
        required=True,
        type=Path,
        help="Output directory for figures.",
    )

    parser.add_argument(
        "--tree",
        type=Path,
        default=None,
        help="Optional Newick tree to order genomes on x-axis.",
    )

    parser.add_argument(
        "--methods",
        nargs="*",
        default=None,
        help="Optional list of methods to plot.",
    )

    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    df_gene_number = read_table(args.tables_dir, "gene_number.tsv")
    df_protein_lengths = read_table(args.tables_dir, "protein_lengths.tsv")
    df_busco = read_table(args.tables_dir, "busco.tsv")
    df_omark = read_table(args.tables_dir, "omark.tsv")
    df_gff = read_table(args.tables_dir, "gff_stats.tsv")

    all_genomes = sorted(set().union(
        set(df_gene_number.get("genome_id", [])),
        set(df_protein_lengths.get("genome_id", [])),
        set(df_busco.get("genome_id", [])),
        set(df_omark.get("genome_id", [])),
        set(df_gff.get("genome_id", [])) if not df_gff.empty else set(),
    ))

    genome_order = get_genome_order(all_genomes, args.tree)
    print(f"Genome order: {genome_order}")

    if args.methods:
        methods = args.methods
    else:
        available_methods = get_methods(
            df_gene_number,
            df_protein_lengths,
            df_busco,
            df_omark,
            df_gff,
        )

        methods = [
            method for method in METHOD_ORDER
            if method in available_methods
        ]

        extra_methods = [
            method for method in available_methods
            if method not in METHOD_ORDER
        ]

        methods.extend(sorted(extra_methods))

    available_pairs = set()

    for df in [
        df_gene_number,
        df_protein_lengths,
        df_busco,
        df_omark,
        df_gff,
    ]:
        if not df.empty and {"genome_id", "method_annotation"}.issubset(df.columns):
            available_pairs.update(
                zip(
                    df["genome_id"].astype(str),
                    df["method_annotation"].astype(str),
                )
            )

    make_combined_figure(
        genome_order,
        methods,
        available_pairs,
        args.out_dir,
        df_gene_number,
        df_protein_lengths,
        df_busco,
        df_omark,
        df_gff,
    )


if __name__ == "__main__":
    main()