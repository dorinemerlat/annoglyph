#!/usr/bin/env python3

from __future__ import annotations
from annotation_phylogeny import load_phylogeny

import argparse
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from Bio import Phylo
from scipy.stats import pearsonr, spearmanr


NEUTRAL_BAR_COLOR = "#8c8c8c"
NEUTRAL_BOX_COLOR = "#bdbdbd"

MAIN_TITLE_SIZE = 13
LABEL_SIZE = 8
TITLE_SIZE = 10
LEGEND_SIZE = 8
LEGEND_TITLE_SIZE = 8

FUNCTIONAL_COLORS = {
    "product_only": "#17C3B2",
    "go_only": "#f0e442",
    "product_and_go": "#69c4ec",
}

plt.rcParams["hatch.linewidth"] = 0.4


def clean_genome_name(genome_id: str) -> str:
    name = re.sub(r"-\d$", "", str(genome_id))
    name = name.replace("-", " ")
    parts = name.split()
    if parts:
        parts[0] = parts[0].capitalize()
    return " ".join(parts)


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Table not found: {path}")

    if path.name.endswith(".stats"):
        return pd.read_csv(path, sep=r"\s+", engine="python")

    return pd.read_csv(path, sep="\t")


def clean_group_name(group_name: str) -> str:
    return str(group_name).split("|")[0]


def get_genome_order(all_genomes, tree_file, taxonomy_file=None):
    if tree_file is None:
        return sorted(all_genomes)
    _, _, order, missing = load_phylogeny(tree_file, all_genomes, taxonomy_file)
    print(f"Found {len(all_genomes) - len(missing)} genomes in tree")
    print(f"Missing from tree: {missing}")
    return order


def get_tree_groups(tree_file: Path, available_genomes: set[str], min_size: int = 4, taxonomy_file=None) -> dict[str, list[str]]:
    if tree_file is None:
        return {}

    tree, leaf_genomes, _, _ = load_phylogeny(tree_file, available_genomes, taxonomy_file)
    groups = {}
    seen = set()

    for i, clade in enumerate(tree.find_clades(order="level")):
        leaves = [genome for leaf in clade.get_terminals()
                  for genome in leaf_genomes[id(leaf)]]

        leaves = sorted(set(leaves))

        if len(leaves) < min_size:
            continue

        if len(leaves) == len(available_genomes):
            group_name = "all"
        else:
            group_name = clean_group_name(clade.name) if clade.name else f"group_{i:03d}"

        key = tuple(leaves)
        if key in seen:
            continue

        seen.add(key)
        groups[group_name] = leaves

    return groups


def style_axis(ax) -> None:
    ax.tick_params(axis="both", labelsize=LABEL_SIZE)
    ax.xaxis.label.set_size(LABEL_SIZE)
    ax.yaxis.label.set_size(LABEL_SIZE)
    ax.xaxis.label.set_weight("bold")
    ax.yaxis.label.set_weight("bold")
    ax.title.set_size(TITLE_SIZE)
    ax.title.set_weight("bold")


def style_legend(legend) -> None:
    if legend is None:
        return

    legend.get_title().set_fontsize(LEGEND_TITLE_SIZE)
    legend.get_title().set_fontweight("bold")

    for text in legend.get_texts():
        text.set_fontsize(LEGEND_SIZE)


def setup_genome_xaxis(ax, genome_order: list[str]) -> None:
    ax.set_xticks(range(len(genome_order)))
    ax.set_xticklabels(
        [clean_genome_name(g) for g in genome_order],
        rotation=90,
        fontsize=LABEL_SIZE,
        fontstyle="italic",
    )


def select_annotation_rows(summary, introns, selection):
    keys = ['genome_id', 'method_annotation']
    selected = selection[keys]
    if selected['genome_id'].duplicated().any():
        raise ValueError('Expected one selected annotation per genome')
    chosen = selected.merge(summary, on=keys, how='left', validate='one_to_one', indicator=True)
    if (chosen['_merge'] != 'both').any():
        raise ValueError('Selected annotations are missing from the summary')
    return chosen.drop(columns='_merge'), introns.merge(selected, on=keys, validate='many_to_one')


def prepare_summary(summary: pd.DataFrame, genome_stats: pd.DataFrame) -> pd.DataFrame:
    df = summary.copy()

    genome_stats = genome_stats.copy()
    genome_stats = genome_stats.rename(columns={
        "sum_len": "assembly_size_bp",
        "N50": "assembly_n50",
        "num_seqs": "assembly_num_seqs",
    })

    keep_cols = [
        "genome_id",
        "assembly_size_bp",
        "assembly_n50",
        "assembly_num_seqs",
    ]

    genome_stats = genome_stats[[c for c in keep_cols if c in genome_stats.columns]]

    df = df.merge(genome_stats, on="genome_id", how="left")

    df["assembly_size_gb"] = df["assembly_size_bp"] / 1e9
    df["assembly_size_kb"] = df["assembly_size_bp"] / 1e3
    df["gene_density_per_kb"] = df["gene_number"] / df["assembly_size_kb"]

    return df


def prepare_intron_stats(intron_lengths: pd.DataFrame) -> pd.DataFrame:
    if intron_lengths.empty:
        return pd.DataFrame(columns=[
            "genome_id",
            "mean_intron_length",
            "median_intron_length",
            "intron_number",
        ])

    return (
        intron_lengths
        .groupby("genome_id", as_index=False)
        .agg(
            mean_intron_length=("intron_length", "mean"),
            median_intron_length=("intron_length", "median"),
            intron_number=("intron_length", "count"),
        )
    )


def plot_summary_figure(
    df: pd.DataFrame,
    intron_lengths: pd.DataFrame,
    genome_order: list[str],
    out_prefix: Path,
) -> None:
    df = df.set_index("genome_id").reindex(genome_order).reset_index()

    n_genomes = len(genome_order)
    fig_width = max(12, n_genomes * 0.38)

    fig, axes = plt.subplots(
        nrows=4,
        ncols=1,
        figsize=(fig_width, 11.69),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 1.2, 1, 1.2]},
    )

    x = range(len(genome_order))

    axes[0].bar(
        x,
        df["mean_exons_per_gene"],
        color=NEUTRAL_BAR_COLOR,
        edgecolor="none",
        width=0.8,
    )
    axes[0].set_title("Mean number of exons per gene", fontweight="bold")
    axes[0].set_ylabel("Mean exons/gene")

    box_data = []
    box_positions = []

    for i, genome_id in enumerate(genome_order):
        values = intron_lengths.loc[
            intron_lengths["genome_id"] == genome_id,
            "intron_length"
        ].dropna()

        if len(values) > 0:
            box_data.append(values)
            box_positions.append(i)

    if box_data:
        bp = axes[1].boxplot(
            box_data,
            positions=box_positions,
            widths=0.6,
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

    axes[1].set_title("Intron length distribution", fontweight="bold")
    axes[1].set_ylabel("Intron length (bp)")
    axes[1].yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, p: f"{int(v):,}"))

    axes[2].bar(
        x,
        df["gene_density_per_kb"],
        color=NEUTRAL_BAR_COLOR,
        edgecolor="none",
        width=0.8,
    )
    axes[2].set_title("Gene density", fontweight="bold")
    axes[2].set_ylabel("Genes per kb")

    functional_cols = [
        "percent_genes_product_only",
        "percent_genes_go_only",
        "percent_genes_product_and_go",
    ]

    functional_labels = [
        "Description only",
        "GO annotation only",
        "Description + GO annotation",
    ]

    functional_colors = [
        FUNCTIONAL_COLORS["product_only"],
        FUNCTIONAL_COLORS["go_only"],
        FUNCTIONAL_COLORS["product_and_go"],
    ]

    bottom = pd.Series([0.0] * len(df))

    for col, label, color in zip(functional_cols, functional_labels, functional_colors):
        values = df[col].fillna(0)
        axes[3].bar(
            x,
            values,
            bottom=bottom,
            color=color,
            edgecolor="none",
            width=0.8,
            label=label,
        )
        bottom += values

    axes[3].set_title("Functional annotation", fontweight="bold")
    axes[3].set_ylabel("Percentage (%)")
    axes[3].set_ylim(0, 100)
    axes[3].yaxis.set_major_formatter(mticker.PercentFormatter())

    handles, labels = axes[3].get_legend_handles_labels()

    legend = axes[3].legend(
        handles[::-1],
        labels[::-1],
        title="Annotation category",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=LEGEND_SIZE,
        title_fontsize=LEGEND_TITLE_SIZE,
    )
    style_legend(legend)

    for ax in axes:
        style_axis(ax)
        ax.set_xlabel("")
        ax.yaxis.set_label_coords(-0.065, 0.5)

    setup_genome_xaxis(axes[-1], genome_order)
    axes[-1].set_xlabel("Genome", fontsize=LABEL_SIZE, fontweight="bold")

    plt.tight_layout(rect=[0.07, 0.04, 0.98, 0.99])

    out_png = Path(str(out_prefix) + ".png")
    out_pdf = Path(str(out_prefix) + ".pdf")
    out_svg = Path(str(out_prefix) + ".svg")

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {out_png}")
    print(f"Saved {out_pdf}")
    print(f"Saved {out_svg}")


def add_correlation_box(ax, x, y, loc: str = "lower_right") -> None:
    data = pd.DataFrame({"x": x, "y": y}).dropna()

    if len(data) < 3:
        text = "r = NA, p = NA\nρ = NA, p = NA"
    else:
        r, p_r = pearsonr(data["x"], data["y"])
        rho, p_s = spearmanr(data["x"], data["y"])

        text = (
            f"$r$ = {r:.3f}, $p$ = {p_r:.2g}\n"
            f"$\\rho$ = {rho:.3f}, $p$ = {p_s:.2g}"
        )

    if loc == "upper_right":
        x_pos, y_pos, va = 0.97, 0.95, "top"
    else:
        x_pos, y_pos, va = 0.97, 0.05, "bottom"

    ax.text(
        x_pos,
        y_pos,
        text,
        transform=ax.transAxes,
        ha="right",
        va=va,
        fontsize=9,
        bbox=dict(
            boxstyle="square,pad=0.3",
            facecolor="white",
            edgecolor="black",
            linewidth=0.8,
        ),
    )


def plot_correlation_panel(
    ax,
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    xlabel: str,
    ylabel: str,
    title: str,
    box_loc: str = "lower_right",
) -> None:
    plot_df = df[[x_col, y_col]].dropna()

    if len(plot_df) == 0:
        ax.set_title(title, fontweight="bold")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        return

    sns.regplot(
        data=plot_df,
        x=x_col,
        y=y_col,
        ax=ax,
        scatter_kws={"s": 18, "color": "black"},
        line_kws={"color": "#2f6fff", "linewidth": 2},
        ci=95,
        color="#2f6fff",
    )

    ax.set_title(title, fontweight="bold")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, p: f"{int(v):,}"))

    add_correlation_box(ax, plot_df[x_col], plot_df[y_col], loc=box_loc)
    style_axis(ax)


def plot_correlation_figure(df: pd.DataFrame, group_name: str, out_prefix: Path) -> None:
    fig, axes = plt.subplots(
        nrows=1,
        ncols=4,
        figsize=(16, 4),
    )

    plot_correlation_panel(
        axes[0],
        df,
        x_col="mean_intron_length",
        y_col="assembly_size_gb",
        xlabel="Mean intron length (bp)",
        ylabel="Assembly size (Gb)",
        title="A",
        box_loc="lower_right",
    )

    plot_correlation_panel(
        axes[1],
        df,
        x_col="intron_number",
        y_col="assembly_size_gb",
        xlabel="Number of introns",
        ylabel="Assembly size (Gb)",
        title="B",
        box_loc="upper_right",
    )

    plot_correlation_panel(
        axes[2],
        df,
        x_col="mean_introns_per_gene",
        y_col="assembly_size_gb",
        xlabel="Mean introns per gene",
        ylabel="Assembly size (Gb)",
        title="C",
        box_loc="upper_right",
    )

    plot_correlation_panel(
        axes[3],
        df,
        x_col="gene_number",
        y_col="assembly_size_gb",
        xlabel="Number of genes",
        ylabel="Assembly size (Gb)",
        title="D",
        box_loc="upper_right",
    )

    fig.suptitle(group_name, fontsize=MAIN_TITLE_SIZE, fontweight="bold")
    plt.tight_layout()

    out_png = Path(str(out_prefix) + ".png")
    out_pdf = Path(str(out_prefix) + ".pdf")
    out_svg = Path(str(out_prefix) + ".svg")

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {out_png}")
    print(f"Saved {out_pdf}")
    print(f"Saved {out_svg}")


def safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate gene-structure summary and correlation figures."
    )

    parser.add_argument(
        "-s",
        "--summary",
        required=True,
        type=Path,
        help="Annotation summary TSV file.",
    )

    parser.add_argument(
        "-i",
        "--intron-lengths",
        required=True,
        type=Path,
        help="Intron lengths TSV file.",
    )

    parser.add_argument(
        "-g",
        "--genome-stats",
        required=True,
        type=Path,
        help="Genome assembly statistics table.",
    )

    parser.add_argument(
        "-o",
        "--out-dir",
        required=True,
        type=Path,
        help="Output directory.",
    )

    parser.add_argument(
        "--tree",
        type=Path,
        default=None,
        help="Optional Newick tree.",
    )

    parser.add_argument(
        "--min-group-size",
        type=int,
        default=4,
        help="Minimum number of genomes per tree group for correlation plots.",
    )

    parser.add_argument("--taxonomy", type=Path, help="Taxonomy TSV mapping scientific names to genome IDs.")
    parser.add_argument('--selection', type=Path, help='Selected annotation method per genome.')
    args = parser.parse_args()

    summary_dir = args.out_dir / "summary_plots"
    correlation_dir = args.out_dir / "correlation_plots"

    summary_dir.mkdir(parents=True, exist_ok=True)
    correlation_dir.mkdir(parents=True, exist_ok=True)

    df_summary = read_table(args.summary)
    df_introns = read_table(args.intron_lengths)
    df_genome_stats = read_table(args.genome_stats)

    if args.selection:
        df_summary, df_introns = select_annotation_rows(df_summary, df_introns, read_table(args.selection))

    intron_stats = prepare_intron_stats(df_introns)

    df = prepare_summary(df_summary, df_genome_stats)
    df = df.merge(intron_stats, on="genome_id", how="left")

    df["mean_introns_per_gene"] = df["intron_number"] / df["gene_number"]

    all_genomes = sorted(set(df["genome_id"].dropna().astype(str)))
    genome_order = get_genome_order(all_genomes, args.tree, args.taxonomy)

    plot_summary_figure(
        df=df,
        intron_lengths=df_introns,
        genome_order=genome_order,
        out_prefix=summary_dir / "gene_structure_summary",
    )

    groups = {"all": genome_order}
    groups.update(get_tree_groups(
        args.tree,
        available_genomes=set(all_genomes),
        min_size=args.min_group_size,
        taxonomy_file=args.taxonomy,
    ))

    for group_name, genomes in groups.items():
        group_df = df[df["genome_id"].isin(genomes)].copy()

        if len(group_df) < args.min_group_size:
            continue

        plot_correlation_figure(
            df=group_df,
            group_name=group_name,
            out_prefix=correlation_dir / f"{safe_name(group_name)}_gene_structure_assembly_correlations",
        )

    if not any(correlation_dir.glob('*.png')):
        # Small test datasets cannot support correlations, but still need a report.
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.axis('off')
        ax.text(0.5, 0.5,
                f'Correlation plots require at least {args.min_group_size} genomes per group.',
                ha='center', va='center', wrap=True)
        for extension in ('png', 'pdf', 'svg'):
            fig.savefig(correlation_dir / f'insufficient_genomes.{extension}', bbox_inches='tight')
        plt.close(fig)


if __name__ == "__main__":
    main()
