#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path
import pandas as pd
import numpy as np

DEFAULT_THRESHOLDS = {'expected_gene_number': 15000, 'gene_number_good_deviation': 5000, 'gene_number_neutral_deviation': 10000, 'protein_length_q1_min': 120, 'busco_complete_min': 70, 'busco_duplicated_max': 3, 'omark_consistent_min': 60, 'omark_contamination_max': 0, 'omark_unknown_max': 30, 'start_codon_min': 80, 'orthoinspector_protein_length_q1_min': 115, 'orthoinspector_busco_complete_min': 70, 'orthofinder_protein_length_q1_min': 150, 'orthofinder_busco_complete_min': 80, 'orthofinder_gene_number_min': 10000, 'orthofinder_contamination_max': 0, 'orthofinder_start_codon_min': 40}

def calculate_final_score(df):
    score_columns = [col for col in df.columns if col.startswith("score_")]
    df["final_score"] = df[score_columns].sum(axis=1)

    return df

def add_scores(df, thresholds):
    target = thresholds['expected_gene_number']
    df['diff_gene_number'] = abs(df['gene_number'] - target) if target is not None else 0
    df['score_gene_number'] = df['diff_gene_number'].apply(
        lambda distance: 1 if distance <= thresholds['gene_number_good_deviation'] else (
            0 if distance < thresholds['gene_number_neutral_deviation'] else -1)
    ) if target is not None else 0
    comparisons = {
        'protein_length_q1': ('q1_protein_length', 'protein_length_q1_min', True),
        'busco_complete': ('busco_complete_percent', 'busco_complete_min', True),
        'busco_duplicated': ('busco_duplicated_percent', 'busco_duplicated_max', False),
        'omark_consistent': ('omark_consistent_percent', 'omark_consistent_min', True),
        'omark_contamination': ('omark_contamination_percent', 'omark_contamination_max', False),
        'omark_unknown': ('omark_unknown_percent', 'omark_unknown_max', False),
    }
    for criterion, (column, key, minimum) in comparisons.items():
        passed = df[column] >= thresholds[key] if minimum else df[column] <= thresholds[key]
        df['score_' + criterion] = passed.astype(int)
    df['score_start_codon'] = (df['percent_genes_start_only'] + df['percent_genes_start_and_stop'] >= thresholds['start_codon_min']).astype(int)
    return calculate_final_score(df)


def select_best_annotation(df):
    max_scores = df.groupby("genome_id")["final_score"].transform("max")
    best_annotations = df[df["final_score"] == max_scores].copy()
    return best_annotations

def add_ranking_scores(df, thresholds):
    new_df = pd.DataFrame()

    duplicated_df = df[df["genome_id"].duplicated(keep=False)].copy()
    unique_df = df[~df["genome_id"].duplicated(keep=False)].copy()

    new_df = pd.concat([new_df, unique_df], ignore_index=True)

    for genome_id in duplicated_df["genome_id"].unique():
        subdf = duplicated_df[duplicated_df["genome_id"] == genome_id].copy()

        subdf["score_gene_number"] = subdf["diff_gene_number"].rank(ascending=False, method="min") if thresholds["expected_gene_number"] is not None else 0
        subdf["score_protein_length_q1"] = subdf["q1_protein_length"].rank(ascending=True, method="min")
        subdf["score_busco_complete"] = subdf["busco_complete_percent"].rank(ascending=True, method="min")
        subdf["score_busco_duplicated"] = subdf["busco_duplicated_percent"].rank(ascending=False, method="min")
        subdf["score_omark_consistent"] = subdf["omark_consistent_percent"].rank(ascending=True, method="min")
        subdf["score_omark_contamination"] = subdf["omark_contamination_percent"].rank(ascending=False, method="min")
        subdf["score_omark_unknown"] = subdf["omark_unknown_percent"].rank(ascending=False, method="min")
        subdf["score_start_codon"] = (
            subdf["percent_genes_start_only"]
            .add(subdf["percent_genes_start_and_stop"])
            .rank(ascending=True, method="min")
        )

        score_columns = [col for col in subdf.columns if col.startswith("score_")]
        subdf["final_score"] = subdf[score_columns].sum(axis=1)

        new_df = pd.concat([new_df, subdf], ignore_index=True)

    return new_df

def add_busco_complete_as_final_score(df):
    df = df.copy()
    df["final_score"] = df["busco_complete_percent"]
    return df

def main():
    parser = argparse.ArgumentParser(
        description="Select the best annotation per genome from annotation summary metrics."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        help="Input annotation summary TSV.",
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output TSV with best annotation per genome.",
    )

    parser.add_argument('--thresholds', type=Path, help='JSON overrides of selection thresholds.')
    parser.add_argument('--scores-output', type=Path, default=Path('selection_scores.tsv'))
    parser.add_argument('--thresholds-output', type=Path, default=Path('selection_thresholds.json'))
    args = parser.parse_args()
    supplied = json.loads(args.thresholds.read_text()) if args.thresholds else {}
    unknown = set(supplied) - set(DEFAULT_THRESHOLDS)
    if unknown:
        parser.error(f'Unknown thresholds: {sorted(unknown)}')
    thresholds = dict(DEFAULT_THRESHOLDS, **supplied)
    for key, value in thresholds.items():
        if key == 'expected_gene_number' and value is None:
            continue
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value) or value < 0:
            parser.error(f'{key} must be a finite nonnegative number')
        if key.endswith(('_min', '_max')) and 'protein_length' not in key and 'gene_number' not in key and value > 100:
            parser.error(f'{key} must be a percentage between 0 and 100')
    if thresholds['gene_number_neutral_deviation'] <= thresholds['gene_number_good_deviation']:
        parser.error('gene_number_neutral_deviation must exceed gene_number_good_deviation')
    args.thresholds_output.write_text(json.dumps(thresholds, indent=2) + '\n')

    df = pd.read_csv(args.input, sep="\t")
    # ne garde que les colonnes d'interet
    interest_cols = [
        "genome_id",
        "method_annotation",
        "gene_number",
        "mean_protein_length",
        "q1_protein_length",
        "busco_complete_percent",
        "busco_duplicated_percent",
        "omark_consistent_percent",
        "omark_contamination_percent",
        "omark_unknown_percent",
        "percent_genes_start_only",
        "percent_genes_start_and_stop",
    ]
    missing = set(interest_cols) - set(df.columns)
    if missing:
        raise ValueError(f'Missing selection metrics: {sorted(missing)}')
    df = df[interest_cols]
    numeric_columns = [column for column in interest_cols if column not in ('genome_id', 'method_annotation')]
    if df.empty or df[numeric_columns].isna().any().any():
        raise ValueError('Automatic selection requires nonempty, complete QC metrics')

    df = add_scores(df, thresholds)
    score_details = df.copy()

    best_annotations = select_best_annotation(df)

    # si au moins un genome_id duplique dans le tableau, met un if
    if best_annotations["genome_id"].duplicated().any():
        print("Warning: multiple best annotations found for some genome_id. Adding ranking scores to break ties.")
        best_annotations = add_ranking_scores(best_annotations, thresholds)
        best_annotations = select_best_annotation(best_annotations)

    if best_annotations["genome_id"].duplicated().any():
        print("Warning: multiple best annotations still found for some genome_id after adding ranking scores. Keeping all tied annotations.")
        # fonction qui Mettre final_score == busco-complete
        best_annotations = add_busco_complete_as_final_score(best_annotations)
        best_annotations = select_best_annotation(best_annotations)

    # Guarantee exactly one final annotation per genome, even after perfect ties.
    # Existing biological score/rank/BUSCO rules above remain unchanged.
    best_annotations = best_annotations.sort_values(
        ["genome_id", "method_annotation"], kind="stable"
    ).drop_duplicates("genome_id", keep="first").copy()

    # maintenant on va ajouter une colonne sufficient quality: yes or no, no siq1 < 120 ou  busco_complete < 70
    best_annotations["sufficient_quality_for_orthoinspector"] = np.where(
        (best_annotations["q1_protein_length"] < thresholds["orthoinspector_protein_length_q1_min"]) | (best_annotations["busco_complete_percent"] < thresholds["orthoinspector_busco_complete_min"]),
        "no",
        "yes"
    )

    # maintenant on va ajouter une colonne sufficient quality: yes or no, no busco-complete < 80 ou q1_protein_length < 150 ou small, contamination > 0 ou start codon < 40
    best_annotations["sufficient_quality_for_orthofinder"] = np.where(
        (best_annotations["q1_protein_length"] < thresholds["orthofinder_protein_length_q1_min"]) | (best_annotations["busco_complete_percent"] < thresholds["orthofinder_busco_complete_min"]) |
        (best_annotations["gene_number"] < thresholds["orthofinder_gene_number_min"]) | (best_annotations["omark_contamination_percent"].fillna(0) > thresholds["orthofinder_contamination_max"]) |
        (best_annotations["percent_genes_start_only"].add(best_annotations["percent_genes_start_and_stop"]) < thresholds["orthofinder_start_codon_min"]),
        "no",
        "yes"
    )
    chosen = set(zip(best_annotations.genome_id, best_annotations.method_annotation))
    score_details['selected'] = [(genome, method) in chosen for genome, method in zip(score_details.genome_id, score_details.method_annotation)]
    score_details.to_csv(args.scores_output, sep='\t', index=False)
    # ne garder dans le tableau que les colonnes genome_id, method_annotation, sufficient_quality_for_orthoinspector, sufficient_quality_for_orthofinder
    best_annotations = best_annotations[["genome_id", "method_annotation", "sufficient_quality_for_orthoinspector", "sufficient_quality_for_orthofinder"]]

    best_annotations.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
