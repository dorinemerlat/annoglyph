#!/usr/bin/env python3

import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Phylo
from matplotlib.colors import LogNorm


AA_ORDER = [
    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Summarize and plot tRNAscan-SE tabular outputs across species.")
    parser.add_argument("--input-dir", type=Path, required=True, help="Root directory containing one subdirectory per species.")
    parser.add_argument("--tree", type=Path, required=True, help="Newick tree containing Diplopoda and Chilopoda clades.")
    parser.add_argument("--seqkit-stats", type=Path, required=True, help="Combined seqkit stats file containing genome_id and sum_len.")
    parser.add_argument("--output-dir", type=Path, default=Path("trnascan_plots"), help="Output directory (default: trnascan_plots).")
    parser.add_argument("--include-pseudogenes", action="store_true", help="Include predictions marked as pseudo in the isotype count table and heatmap.")
    return parser.parse_args()


def normalize_name(value):
    return re.sub(r"[^a-z0-9]+", "-", str(value).strip().lower()).strip("-")


def canonical_species_id(value):
    species_id = normalize_name(str(value).split("|", 1)[0])
    return species_id[:-2] if species_id.endswith("-1") else species_id


def display_name(species_id):
    return species_id.replace("-", " ").capitalize()


def read_classes_from_tree(path):
    tree = Phylo.read(path, "newick")
    tree_order = list(dict.fromkeys(canonical_species_id(terminal.name) for terminal in tree.get_terminals()))
    tree_species = set(tree_order)
    classes = {}
    for class_name in ["Diplopoda", "Chilopoda"]:
        class_clades = [clade for clade in tree.find_clades() if clade.name and clade.name.split("|", 1)[0] == class_name]
        if len(class_clades) != 1:
            raise ValueError(f"Expected one {class_name} clade in {path}, found {len(class_clades)}.")
        for terminal in class_clades[0].get_terminals():
            classes[canonical_species_id(terminal.name)] = class_name
    return classes, tree_species, tree_order


def read_genome_sizes(path):
    table = pd.read_csv(path, sep=r"\s+")
    required_columns = {"genome_id", "sum_len"}
    if not required_columns.issubset(table.columns):
        raise ValueError(f"{path} must contain the columns genome_id and sum_len.")
    table["species_id"] = table["genome_id"].map(canonical_species_id)
    if table["species_id"].duplicated().any():
        duplicates = ", ".join(sorted(table.loc[table["species_id"].duplicated(False), "species_id"].unique()))
        raise ValueError(f"Duplicated genome identifiers after normalization: {duplicates}")
    return table.set_index("species_id")["sum_len"].astype(int).to_dict()


def select_output_files(input_dir):
    candidates = sorted(input_dir.rglob("*_trna.out"))
    if not candidates:
        raise FileNotFoundError(f"No *_trna.out files found under {input_dir}")

    by_species = {}
    for path in candidates:
        species_id = path.parent.name
        by_species.setdefault(species_id, []).append(path)

    selected = {}
    for species_id, paths in by_species.items():
        expected = [path for path in paths if path.name == f"{species_id}_trna.out"]
        if len(expected) == 1:
            selected[species_id] = expected[0]
        elif len(paths) == 1:
            selected[species_id] = paths[0]
        else:
            choices = "\n  ".join(str(path) for path in paths)
            raise ValueError(f"Several ambiguous tRNAscan-SE outputs were found for {species_id}:\n  {choices}")
    return selected


def select_one_output_per_species(selected_files):
    by_species = {}
    for species_id, path in selected_files.items():
        by_species.setdefault(canonical_species_id(species_id), []).append((species_id, path))

    unique_files = {}
    for species_key, entries in by_species.items():
        if len(entries) == 1:
            species_id, path = entries[0]
            unique_files[species_id] = path
            continue

        nonempty = [(species_id, path) for species_id, path in entries if path.stat().st_size > 0]
        candidates = nonempty if nonempty else entries
        canonical = [(species_id, path) for species_id, path in candidates if normalize_name(species_id) == species_key]

        if len(canonical) == 1:
            chosen = canonical[0]
        elif len(candidates) == 1:
            chosen = candidates[0]
        else:
            choices = "\n  ".join(str(path) for _, path in candidates)
            raise ValueError(f"Several non-empty tRNAscan-SE outputs match {species_key}:\n  {choices}")

        ignored = [str(path) for species_id, path in entries if (species_id, path) != chosen]
        print(f"Warning: duplicate outputs ignored for {species_key}: {', '.join(ignored)}", file=sys.stderr)
        unique_files[chosen[0]] = chosen[1]

    return unique_files


def parse_trnascan_output(path, species_id, species_label, taxonomic_class):
    if path.stat().st_size == 0:
        print(f"Warning: empty tRNAscan-SE output skipped for {species_id}: {path}", file=sys.stderr)
        return []

    records = []
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 9 or not fields[1].isdigit():
                continue
            note = " ".join(fields[9:])
            records.append({
                "species_id": species_id,
                "species": species_label,
                "class": taxonomic_class,
                "sequence": fields[0],
                "trna_number": int(fields[1]),
                "begin": int(fields[2]),
                "end": int(fields[3]),
                "isotype": fields[4],
                "anticodon": fields[5],
                "intron_begin": int(fields[6]),
                "intron_end": int(fields[7]),
                "score": float(fields[8]),
                "note": note,
                "is_pseudogene": "pseudo" in note.lower(),
                "source": str(path),
            })
    if not records:
        print(f"Warning: no valid tRNA prediction found for {species_id}; file skipped: {path}", file=sys.stderr)
    return records


def save_figure(fig, output_stem):
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(counts, species_order, output_dir):
    matrix = counts.pivot(index="species", columns="isotype", values="count").fillna(0)
    matrix = matrix.reindex(columns=AA_ORDER, fill_value=0)
    ordered_labels = [display_name(species_id) for species_id in species_order]
    matrix = matrix.reindex([label for label in ordered_labels if label in matrix.index])

    fig_height = max(4.5, 0.30 * len(matrix))
    fig, ax = plt.subplots(figsize=(12, fig_height))
    positive = matrix.to_numpy()[matrix.to_numpy() > 0]
    use_log_scale = positive.size > 0 and positive.max() / positive.min() >= 100
    norm = LogNorm(vmin=max(1, positive.min()), vmax=positive.max()) if use_log_scale else None
    sns.heatmap(matrix, cmap="viridis", norm=norm, linewidths=0.15, linecolor="white", cbar_kws={"label": "tRNA count"}, ax=ax)
    ax.set_xlabel("Amino acid isotype")
    ax.set_ylabel("")
    ax.set_title("Predicted tRNA genes by amino acid isotype")
    fig.tight_layout()
    save_figure(fig, output_dir / "trna_isotype_heatmap")


def plot_totals(summary, species_order, output_dir):
    order = {species_id: index for index, species_id in enumerate(species_order)}
    data = summary.copy()
    data["tree_order"] = data["species_id"].map(lambda value: order[canonical_species_id(value)])
    data = data.sort_values("tree_order", ascending=False)
    fig_height = max(4.5, 0.30 * len(data))
    fig, ax = plt.subplots(figsize=(10, fig_height))
    ax.barh(data["species"], data["functional_trnas"], label="Functional predictions", color="#0072B2")
    ax.barh(data["species"], data["pseudogenes"], left=data["functional_trnas"], label="Pseudogenes", color="#D55E00")
    ax.set_xlabel("Number of predicted tRNAs")
    ax.set_ylabel("")
    ax.set_title("tRNAscan-SE predictions per species")
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    save_figure(fig, output_dir / "trna_total_counts")


def plot_violin(summary, output_dir):
    data = summary[summary["class"].isin(["Chilopoda", "Diplopoda"])].copy()
    class_counts = data.groupby("class")["species"].nunique()
    if len(class_counts) < 2 or (class_counts < 2).any():
        print("Warning: at least two species per class are required for the violin plot; plot skipped.", file=sys.stderr)
        return

    order = [group for group in ["Chilopoda", "Diplopoda"] if group in class_counts.index]
    metrics = [
        ("functional_trnas", "Predicted functional tRNAs", "Counts"),
        ("functional_trnas_per_mb", "Predicted functional tRNAs per Mb", "Density"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
    for ax, (metric, ylabel, title) in zip(axes, metrics):
        metric_data = data.dropna(subset=[metric])
        sns.violinplot(data=metric_data, x="class", y=metric, order=order, inner=None, cut=0, color="#9ECAE1", ax=ax)
        sns.stripplot(data=metric_data, x="class", y=metric, order=order, color="#222222", size=5, jitter=0.12, ax=ax)
        values = metric_data[metric]
        if not values.empty and values.min() > 0 and values.max() / values.min() >= 100:
            ax.set_yscale("log")
            ylabel += " (log scale)"
        ax.set_xlabel("")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        sns.despine(ax=ax)
    fig.suptitle("tRNA gene predictions by myriapod class", y=1.02)
    fig.tight_layout()
    save_figure(fig, output_dir / "trna_counts_by_class")


def plot_genome_size(summary, output_dir):
    data = summary.dropna(subset=["genome_size_mb"]).copy()
    data = data[data["class"].isin(["Chilopoda", "Diplopoda"])]
    if data.empty:
        print("Warning: no matched genome sizes; genome-size plot skipped.", file=sys.stderr)
        return

    fig, ax = plt.subplots(figsize=(7.5, 6))
    sns.scatterplot(data=data, x="genome_size_mb", y="functional_trnas", hue="class", style="class", s=65, ax=ax)
    for _, row in data.iterrows():
        ax.annotate(row["species"], (row["genome_size_mb"], row["functional_trnas"]), xytext=(4, 3), textcoords="offset points", fontsize=6)
    if data["functional_trnas"].min() > 0 and data["functional_trnas"].max() / data["functional_trnas"].min() >= 100:
        ax.set_yscale("log")
        ax.set_ylabel("Predicted functional tRNAs (log scale)")
    else:
        ax.set_ylabel("Predicted functional tRNAs")
    ax.set_xlabel("Assembly size (Mb)")
    ax.set_title("Assembly size and predicted tRNA content")
    ax.legend(title="Class", frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    save_figure(fig, output_dir / "trna_counts_vs_genome_size")


def main():
    args = parse_args()
    classes, tree_species, tree_order = read_classes_from_tree(args.tree)
    genome_sizes = read_genome_sizes(args.seqkit_stats)
    selected_files = select_output_files(args.input_dir)
    outside_tree = sorted(species_id for species_id in selected_files if canonical_species_id(species_id) not in tree_species)
    selected_files = {
        species_id: path
        for species_id, path in selected_files.items()
        if canonical_species_id(species_id) in tree_species
    }
    selected_files = select_one_output_per_species(selected_files)
    if outside_tree:
        print(f"Warning: outputs excluded because they are absent from the tree: {', '.join(outside_tree)}", file=sys.stderr)
    if not selected_files:
        raise ValueError("No tRNAscan-SE output corresponds to a species present in the tree.")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    records = []
    for species_id, path in selected_files.items():
        species_key = canonical_species_id(species_id)
        species_label = display_name(species_key)
        taxonomic_class = classes.get(species_key, "Unknown")
        records.extend(parse_trnascan_output(path, species_id, species_label, taxonomic_class))

    if not records:
        raise ValueError("No valid tRNA predictions were found in the selected output files.")

    predictions = pd.DataFrame(records)
    selected_predictions = predictions if args.include_pseudogenes else predictions[~predictions["is_pseudogene"]]

    counts = selected_predictions.groupby(["species_id", "species", "class", "isotype"], as_index=False).size().rename(columns={"size": "count"})
    summary = predictions.groupby(["species_id", "species", "class"], as_index=False).agg(
        total_predictions=("trna_number", "size"),
        pseudogenes=("is_pseudogene", "sum"),
    )
    summary["functional_trnas"] = summary["total_predictions"] - summary["pseudogenes"]
    summary["genome_size_bp"] = summary["species_id"].map(lambda value: genome_sizes.get(canonical_species_id(value)))
    summary["genome_size_mb"] = summary["genome_size_bp"] / 1_000_000
    summary["functional_trnas_per_mb"] = summary["functional_trnas"] / summary["genome_size_mb"]

    unknown_classes = sorted(summary.loc[summary["class"] == "Unknown", "species"].unique())
    missing_sizes = sorted(summary.loc[summary["genome_size_bp"].isna(), "species"].unique())
    if unknown_classes:
        print(f"Warning: species absent from the Chilopoda/Diplopoda clades: {', '.join(unknown_classes)}", file=sys.stderr)
    if missing_sizes:
        print(f"Warning: species absent from seqkit stats: {', '.join(missing_sizes)}", file=sys.stderr)

    predictions.to_csv(args.output_dir / "trna_predictions.tsv", sep="\t", index=False)
    counts.to_csv(args.output_dir / "trna_counts_by_isotype.tsv", sep="\t", index=False)
    summary.to_csv(args.output_dir / "trna_summary.tsv", sep="\t", index=False)

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.05)
    standard_counts = counts[counts["isotype"].isin(AA_ORDER)]
    plot_heatmap(standard_counts, tree_order, args.output_dir)
    plot_totals(summary, tree_order, args.output_dir)
    plot_violin(summary, args.output_dir)
    plot_genome_size(summary, args.output_dir)

    processed_species = predictions["species_id"].nunique()
    skipped_species = len(selected_files) - processed_species
    print(f"Processed {processed_species} species; skipped {skipped_species} empty or invalid outputs.")
    print(f"Results written to {args.output_dir}")


if __name__ == "__main__":
    main()
