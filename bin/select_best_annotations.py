#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
import numpy as np

def calculate_final_score(df):
    score_columns = [col for col in df.columns if col.startswith("score_")]
    df["final_score"] = df[score_columns].sum(axis=1)

    return df 

def add_scores(df):
    # ajoute les colonnes score_gene_number: 1 si abs (gene_number - 15000) <= 5000, 0 si abs (gene_number - 15000) > 5000 et < 10000, -1 si abs (gene_number - 15000) == 5000
    df["diff_gene_number"] = df["gene_number"].apply(lambda x: abs(x - 15000))
    df["score_gene_number"] = df["diff_gene_number"].apply(lambda x: 1 if x <= 5000 else (0 if x > 5000 and x < 10000 else -1))

    # ajoute la colonne score_protein_length_q1: 1 si q1_protein_length >= 120, sinon 0
    df["score_protein_length_q1"] = df["q1_protein_length"].apply(lambda x: 1 if x >= 120 else 0)

    # ajoute la colonne score_busco_complete: 1 si busco_complete_percent >= 70, sinon 0 
    df["score_busco_complete"] = df["busco_complete_percent"].apply(lambda x: 1 if x >= 70 else 0)

    # ajoute la colonne score_busco_duplicated: 1 si busco_duplicated_percent <= 3, sinon 0
    df["score_busco_duplicated"] = df["busco_duplicated_percent"].apply(lambda x: 1 if x <= 3 else 0)

    # ajoute la colonne score_omark_consistent: 1 si omark_consistent_percent >= 60, sinon 0
    df["score_omark_consistent"] = df["omark_consistent_percent"].apply(lambda x: 1 if x >= 60 else 0)

    # ajoute la colonne score_omark_contamination: 1 si omark_contamination_percent == 0, sinon 0
    df["score_omark_contamination"] = df["omark_contamination_percent"].apply(lambda x: 1 if x == 0 else 0)

    # ajoute la colonne score_omark_unknown: 1 si omark_unknown_percent <= 30, sinon 0
    df["score_omark_unknown"] = df["omark_unknown_percent"].apply(lambda x: 1 if x <= 30 else 0)

    # ajoute la colonne score_start_codon: 1 si percent_genes_start_only + percent_genes_start_and_stop >= 80, sinon 0
    df["score_start_codon"] = df["percent_genes_start_only"].add(df["percent_genes_start_and_stop"]).apply(lambda x: 1 if x >= 80 else 0)

    # final score: somme de toutes les colonnes score_*
    df = calculate_final_score(df)

    return df


def select_best_annotation(df):
    max_scores = df.groupby("genome_id")["final_score"].transform("max")
    best_annotations = df[df["final_score"] == max_scores].copy()
    return best_annotations

def add_ranking_scores(df):
    new_df = pd.DataFrame()

    duplicated_df = df[df["genome_id"].duplicated(keep=False)].copy()
    unique_df = df[~df["genome_id"].duplicated(keep=False)].copy()

    new_df = pd.concat([new_df, unique_df], ignore_index=True)

    for genome_id in duplicated_df["genome_id"].unique():
        subdf = duplicated_df[duplicated_df["genome_id"] == genome_id].copy()

        subdf["score_gene_number"] = subdf["diff_gene_number"].rank(ascending=False, method="min")
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

    args = parser.parse_args()

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
    df = df[interest_cols]

    df = add_scores(df)

    best_annotations = select_best_annotation(df)

    # si au moins un genome_id duplique dans le tableau, met un if 
    if best_annotations["genome_id"].duplicated().any():
        print("Warning: multiple best annotations found for some genome_id. Adding ranking scores to break ties.")
        best_annotations = add_ranking_scores(best_annotations)
        best_annotations = select_best_annotation(best_annotations)
    
    if best_annotations["genome_id"].duplicated().any():
        print("Warning: multiple best annotations still found for some genome_id after adding ranking scores. Keeping all tied annotations.")
        # fonction qui Mettre final_score == busco-complete 
        best_annotations = add_busco_complete_as_final_score(best_annotations)
        best_annotations = select_best_annotation(best_annotations)

    # maintenant on va ajouter une colonne sufficient quality: yes or no, no siq1 < 120 ou  busco_complete < 70
    best_annotations["sufficient_quality_for_orthoinspector"] = np.where(
        (best_annotations["q1_protein_length"] < 115) | (best_annotations["busco_complete_percent"] < 70),
        "no",
        "yes"
    )

    # maintenant on va ajouter une colonne sufficient quality: yes or no, no busco-complete < 80 ou q1_protein_length < 150 ou small, contamination > 0 ou start codon < 40
    best_annotations["sufficient_quality_for_orthofinder"] = np.where(
        (best_annotations["q1_protein_length"] < 150) | (best_annotations["busco_complete_percent"] < 80) |
        (best_annotations["gene_number"] < 10000) | (best_annotations["omark_contamination_percent"].fillna(0) > 0) | 
        (best_annotations["percent_genes_start_only"].add(best_annotations["percent_genes_start_and_stop"]) < 40),
        "no",
        "yes"
    )
    # ne garder dans le tableau que les colonnes genome_id, method_annotation, sufficient_quality_for_orthoinspector, sufficient_quality_for_orthofinder
    best_annotations = best_annotations[["genome_id", "method_annotation", "sufficient_quality_for_orthoinspector", "sufficient_quality_for_orthofinder"]]

    best_annotations.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()