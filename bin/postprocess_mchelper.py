#!/usr/bin/env python3
"""Post-process MCHelper results and hmmscan hits.

This script:
- reads the family table from ``clean_repeat_families.py``,
- adds the best gene and rRNA hmmscan hits,
- merges MCHelper FASTA and .classif information,
- filters out sequences that look like gene/rRNA contaminants,
- and writes a final annotated family table.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    p = argparse.ArgumentParser(
        description="Post-process MCHelper hmmscan outputs and merge into families table"
    )
    p.add_argument('--table', type=Path, required=True, help='Families TSV from clean_repeat_families.py')
    p.add_argument('--fasta', type=Path, required=True, help='FASTA prepared for MCHelper (cleaned_for_mchelper.fa)')
    p.add_argument('--mchelper-classified-fasta', dest='mchelper_classified_fasta', type=Path, required=True, help='MCHelper classified FASTA (classifiedModule/kept_seqs_classified_module.fa)')
    p.add_argument('--mchelper-unclassified-fasta', dest='mchelper_unclassified_fasta', type=Path, required=True, help='MCHelper unclassified FASTA (unclassifiedModule/kept_seqs_unclassified_module.fa)')
    p.add_argument('--mchelper-classif', dest='mchelper_classif', type=Path, required=True, help='MCHelper .classif file for classified module (classifiedModule/denovoLibTEs_PC.classif)')
    p.add_argument('--hmmscan-genes', dest='hmmscan_genes', type=Path, required=True, help='hmmscan tblout-like file for genes')
    p.add_argument('--hmmscan-rrna', dest='hmmscan_rrna', type=Path, required=True, help='hmmscan tblout-like file for rRNA')
    p.add_argument('--out', dest='out_prefix', type=Path, required=True, help='Output prefix: writes <prefix>.processed.tsv')
    p.add_argument('--classification', dest='classification', type=Path, required=True, help='Classification CSV used to replace "Unknown" rm_code values when MCHelper is present (columns: Code, Class, Order, SuperFamily)')
    return p.parse_args()


def validate_inputs(families_path: Path, fasta_path: Path, mchelper_classified_fasta: Path, mchelper_unclassified_fasta: Path, mchelper_classif: Path, genes_path: Path, rrna_path: Path, classification_path: Path) -> None:
    """Fail fast if required inputs are missing."""

    required = [
        (families_path, "families table"),
        (fasta_path, "MCHelper FASTA"),
        (mchelper_classified_fasta, "MCHelper classified FASTA"),
        (mchelper_unclassified_fasta, "MCHelper unclassified FASTA"),
        (mchelper_classif, "MCHelper classif file"),
        (genes_path, "hmmscan genes output"),
        (rrna_path, "hmmscan rRNA output"),
        (classification_path, "classification CSV"),
    ]

    missing: list[str] = []
    for path, label in required:
        if path is None or not path.exists():
            missing.append(f"{label}: {path}")

    if missing:
        print("ERROR: missing required inputs:", file=sys.stderr)
        for item in missing:
            print(f"  - {item}", file=sys.stderr)
        sys.exit(2)


def parse_hmmscan_tbl(path: Path, hit_type: str) -> pd.DataFrame:
    """Parse a hmmscan tblout-like file into a small DataFrame.

    Uses only: target (column 0), query (column 2), e-value (4), score (5).
    Query IDs are trimmed at the first ``#`` to match TE family IDs.
    """
    rows = []
    with path.open('rt') as fh:
        for ln in fh:
            ln = ln.rstrip('\n')
            if not ln or ln.startswith('#'):
                continue
            parts = ln.split()
            if len(parts) < 6:
                continue
            target = parts[0]
            query_full = parts[2]
            query_short = query_full.split('#', 1)[0]
            evalue = parts[4]
            score = parts[5]

            rows.append(
                {
                    "id": query_short,
                    f"{hit_type}_target": target,
                    f"{hit_type}_evalue": evalue,
                    f"{hit_type}_score": score,
                }
            )
    if not rows:
        return pd.DataFrame(
            columns=["id", f"{hit_type}_target", f"{hit_type}_evalue", f"{hit_type}_score"]
        )
    df = pd.DataFrame(rows)

    # Force datatypes
    df = df.astype(
        {
            "id": "string",
            f"{hit_type}_target": "string",
            f"{hit_type}_evalue": "float64",
            f"{hit_type}_score": "float64",
        }
    )

    return df


def parse_fasta_lengths(path: Path) -> pd.DataFrame:
    """Read MCHelper FASTA and extract id, classification, length and comment."""
    rows = []
    records = list(SeqIO.parse(path, "fasta"))
    for record in records:
        base_id, classification = record.id.split("#", 1)

        parts = base_id.split("_")
        if len(parts) > 2:
            # Extract comment after first two components of the ID
            comment = "_".join(parts[2:])
            new_base_id = "_".join(parts[:2])
        else:
            new_base_id = base_id
            comment = ""

        rows.append(
            {
                "id": new_base_id,
                'mchelper_id': base_id,
                "MCHelper": classification,
                "extended_length": len(record.seq),
                "comment": comment,
                "replace_sequence": True,
            }
        )
    df = pd.DataFrame(rows)
    df = df.astype(
        {
            "id": "string",
            'mchelper_id': "string",
            "MCHelper": "string",
            "extended_length": "int64",
            "comment": "string",
            "replace_sequence": "bool",
        }
    )

    return records, df


def parse_classif_file(path: Path) -> pd.DataFrame:
    """Parse a MCHelper ``.classif`` file and normalize identifiers.

    - ``Seq_name`` / ``seq_name`` → ``id``
    - IDs are collapsed from ``A_B_xxx`` to ``A_B`` and the suffix becomes
      a ``comment`` when present.
    - All columns except ``id`` and ``comment`` are prefixed with
      ``MCHelper_`` to avoid clashes with other tables.
    """

    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    # Normalize ID column name
    df = df.rename(columns={"Seq_name": "id", "seq_name": "id"})

    # If ID column is missing, create an empty one and keep going
    if "id" not in df.columns:
        df["id"] = pd.Series(dtype="string")
        df["comment"] = ""
        # don't return here: we still want column prefixing below

    # Normalize IDs and extract comments from suffixes
    new_ids: list[str] = []
    comments: list[str] = []

    for raw_id in df["id"]:
        parts = raw_id.split("_")

        if len(parts) > 2:
            new_id = "_".join(parts[:2])
            comment = "_".join(parts[2:])
        else:
            new_id = raw_id
            comment = ""

        new_ids.append(new_id)
        comments.append(comment)

    df["id"] = pd.Series(new_ids, dtype="string")
    df["comment"] = comments

    # Prefix all column names except 'id' and 'comment' with 'MCHelper_' to avoid collisions
    new_cols: dict[str, str] = {}
    for col in df.columns:
        if col in {"id", "comment"}:
            new_cols[col] = col
        else:
            new_cols[col] = f"MCHelper_{col}"

    df = df.rename(columns=new_cols)

    return df


def parse_classification_csv(path: Path) -> pd.DataFrame:
    """Read a classification CSV and normalize to columns Code, Class, Order, SuperFamily.

    The function is robust to varying column name capitalizations. It returns a
    DataFrame with columns 'Code','Class','Order','SuperFamily' and trims values.
    """
    df = pd.read_csv(path, dtype=str, keep_default_na=False)

    # Normalize column names (case-insensitive lookup)
    cols = {c.lower(): c for c in df.columns}

    def find_col(*names):
        for n in names:
            if n.lower() in cols:
                return cols[n.lower()]
        return None

    code_col = find_col('Code', 'code')
    class_col = find_col('Class', 'class')
    order_col = find_col('Order', 'order')
    super_col = find_col('SuperFamily', 'superfamily', 'Superfamily', 'super_family')

    # Build a new DataFrame with the expected columns (fill missing with empty string)
    new = pd.DataFrame()
    new['Code'] = df[code_col] if code_col is not None else pd.Series([''] * len(df))
    new['Class'] = df[class_col] if class_col is not None else pd.Series([''] * len(df))
    new['Order'] = df[order_col] if order_col is not None else pd.Series([''] * len(df))
    new['SuperFamily'] = df[super_col] if super_col is not None else pd.Series([''] * len(df))

    # Trim whitespace
    for c in ['Code', 'Class', 'Order', 'SuperFamily']:
        new[c] = new[c].astype(str).str.strip()

    # If Code contains multiple values separated by ';', keep the first one
    new['Code'] = new['Code'].astype(str).str.split(';').str[0].str.strip()

    return new


def update_rm_code_from_mchelper(te_table: pd.DataFrame, classification_df: pd.DataFrame) -> pd.DataFrame:
    """Update rows where MCHelper is present AND rm_code == 'Unknown'.

    For each matching row, split MCHelper on '/' to get Class/Order/SuperFamily
    (missing SuperFamily becomes 'Unknown'). Then lookup a matching row in
    classification_df and, if found, replace rm_code with Code and set
    rm2_class/rm2_order/rm2_superfamily with the values from classification_df.

    Matching is case-insensitive and trims whitespace. The function returns a
    modified copy of te_table and prints a short summary to stderr.
    """
    if 'rm_code' not in te_table.columns:
        print("Note: 'rm_code' column not present — skipping classification update", file=sys.stderr)
        return te_table

    # Prepare lookup dict from classification_df: key=(class,order,superfamily) lowercased
    lookup: dict[tuple[str, str, str], dict] = {}
    for _, r in classification_df.iterrows():
        k = (str(r['Class']).strip().lower(), str(r['Order']).strip().lower(), str(r['SuperFamily']).strip().lower())
        lookup[k] = {'Code': r['Code'], 'Class': r['Class'], 'Order': r['Order'], 'SuperFamily': r['SuperFamily']}

    df = te_table.copy()

    # Add 'statut' column: if MCHelper or MCHelper_Reason contains 'MITE' then 'MITE', else empty
    # Be robust to missing columns.
    mch_helper_series = df['MCHelper'] if 'MCHelper' in df.columns else pd.Series([''] * len(df), index=df.index)
    mch_reason_series = df['MCHelper_Reason'] if 'MCHelper_Reason' in df.columns else pd.Series([''] * len(df), index=df.index)
    mite_mask = (
        mch_helper_series.fillna('').astype(str).str.contains('MITE', case=False, na=False)
        | mch_reason_series.fillna('').astype(str).str.contains('MITE', case=False, na=False)
    )
    df['statut'] = ''
    if mite_mask.any():
        df.loc[mite_mask, 'statut'] = 'MITE'

    # If statut == 'MITE' and rm_code is 'Unknown', set conservative defaults
    # for MITE: classify as DNA / ClassII / TIR / Unknown superfamily.
    # Ensure target columns exist first.
    for col in ('rm2_class', 'rm2_order', 'rm2_superfamily'):
        if col not in df.columns:
            df[col] = ''

    mite_default_mask = (df['statut'] == 'MITE') & (df.get('rm_code') == 'Unknown')
    if mite_default_mask.any():
        df.loc[mite_default_mask, 'rm_code'] = 'DNA'
        df.loc[mite_default_mask, 'rm2_class'] = 'ClassII'
        df.loc[mite_default_mask, 'rm2_order'] = 'TIR'
        df.loc[mite_default_mask, 'rm2_superfamily'] = 'Unknown'
        print(f"Applied MITE defaults to {int(mite_default_mask.sum())} rows (rm_code->DNA, rm2_class->ClassII, rm2_order->TIR)", file=sys.stderr)

    mch_present = ~_col_missing_or_empty(df, 'MCHelper')
    rm_unknown = df['rm_code'] == 'Unknown'
    mask = mch_present & rm_unknown

    indices = df.index[mask].tolist()
    if not indices:
        print('No rows to update: either no MCHelper present or no rm_code=="Unknown" rows', file=sys.stderr)
        return df

    updated = 0
    not_found = 0

    # Ensure target columns exist
    for col in ('rm2_class', 'rm2_order', 'rm2_superfamily'):
        if col not in df.columns:
            df[col] = ''

    for idx in indices:
        mch = str(df.at[idx, 'MCHelper']).strip()
        parts = mch.split('/') if mch else []
        cls = parts[0].strip() if len(parts) > 0 else ''
        ordv = parts[1].strip() if len(parts) > 1 else ''
        sup = parts[2].strip() if len(parts) > 2 else 'Unknown'

        key = (cls.lower(), ordv.lower(), sup.lower())
        info = lookup.get(key)
        if info:
            # classification 'Code' may contain multiple values separated by ';'
            code_val = str(info.get('Code', '')).split(';')[0].strip()
            df.at[idx, 'rm_code'] = code_val
            df.at[idx, 'rm2_class'] = info.get('Class', '')
            df.at[idx, 'rm2_order'] = info.get('Order', '')
            df.at[idx, 'rm2_superfamily'] = info.get('SuperFamily', '')
            updated += 1
        else:
            not_found += 1

    print(f"Updated {updated} rows from classification CSV; {not_found} rows had no matching classification", file=sys.stderr)
    return df

def merge_tables(left: pd.DataFrame, right: pd.DataFrame, on: str | list[str] = "id") -> pd.DataFrame:
    """Thin wrapper around :meth:`pandas.DataFrame.merge` using a left join."""

    return left.merge(right, on=on, how="left")


def select_min_evalue(df: pd.DataFrame, hit_type: str) -> pd.DataFrame:
    """For each ``id``, keep the row with the smallest ``<hit_type>_evalue``.

    NaNs are preserved in the final table; they are only replaced by
    ``+inf`` in a temporary column used for the ``idxmin`` computation.
    """

    score_col = f"{hit_type}_evalue"
    tmp_col = f"_tmp_{hit_type}_evalue"
    df_tmp = df.assign(**{tmp_col: df[score_col].fillna(np.inf)})

    best_idx = df_tmp.groupby("id")[tmp_col].idxmin()
    return df.loc[best_idx].reset_index(drop=True)


def _col_missing_or_empty(df: pd.DataFrame, col: str) -> pd.Series:
    """Return boolean Series True when column is missing or value is NA or empty string."""
    if col not in df.columns:
        # If column absent, consider it as missing/empty for all rows
        return pd.Series([True] * len(df), index=df.index)
    s = df[col]
    # treat NaN or empty-string as missing
    return s.isna() | (s.astype(str).str.strip() == "")


def filter_out_gene_rrna_hits(df: pd.DataFrame) -> pd.DataFrame:
    """Filter out rows that are likely gene/rRNA contaminants.

    Removal rule (translated from user's spec):
    - If there is a gene hit (non-empty `gene_target`) OR
      there is a rRNA hit (non-empty `rrna_target`) with a very
      significant e-value (<= 1e-10), and
    - Both `MCHelper_coding` and `MCHelper_struct` are missing/empty,
    then consider this sequence a hit to the gene/rRNA databases and remove it.
    """
    df = df.copy()

    # Detect presence of gene hit
    has_gene = (~_col_missing_or_empty(df, 'gene_target'))

    # Detect significant rrna hit
    has_rrna = (~_col_missing_or_empty(df, 'rrna_target'))
    rrna_e = df['rrna_evalue'] if 'rrna_evalue' in df.columns else pd.Series([np.nan] * len(df), index=df.index)
    rrna_sig = has_rrna & rrna_e.fillna(np.inf).le(1e-10)

    # Detect absent MCHelper annotations
    coding_missing = _col_missing_or_empty(df, 'MCHelper_coding')
    struct_missing = _col_missing_or_empty(df, 'MCHelper_struct')

    remove_mask = (has_gene | rrna_sig) & coding_missing & struct_missing

    df_removed = df.loc[remove_mask].copy()
    n_remove = int(remove_mask.sum())
    if n_remove > 0:
        print(f"Filtering out {n_remove} rows that match gene/rRNA banks (kept {len(df)-n_remove} rows)", file=sys.stderr)
        # Print removed rows as TSV to stderr so user can inspect them
        try:
            df_removed.to_csv(sys.stderr, sep='\t', index=False)
        except Exception:
            # fallback to plain print if writing to stderr as file-like fails
            print(df_removed.to_string(index=False), file=sys.stderr)

    return df.loc[~remove_mask].reset_index(drop=True)


def main() -> None:
    args = parse_args()
    te_file = args.table
    fasta_file = args.fasta
    mchelper_classified_fasta = args.mchelper_classified_fasta
    mchelper_unclassified_fasta = args.mchelper_unclassified_fasta
    mchelper_classif = args.mchelper_classif
    gene_file = args.hmmscan_genes
    rrna_file = args.hmmscan_rrna
    out_prefix = args.out_prefix
    classification = args.classification

    validate_inputs(
        te_file,
        fasta_file,
        mchelper_classified_fasta,
        mchelper_unclassified_fasta,
        mchelper_classif,
        gene_file,
        rrna_file,
        classification,
    )
    te_table_raw = pd.read_csv(te_file, sep='\t', dtype=str)

    # 1. Add informations about genes and rRNA identification from hmmscan outputs
    gene_table = parse_hmmscan_tbl(gene_file, 'gene')
    gene_table_filtered = select_min_evalue(gene_table, 'gene')
    rrna_table = parse_hmmscan_tbl(rrna_file, 'rrna')
    rrna_table_filtered = select_min_evalue(rrna_table, 'rrna')

    # merge genes (LEFT JOIN: keep all te_df rows and add gene info when available)
    te_table = merge_tables(te_table_raw, gene_table_filtered)
    te_table = merge_tables(te_table, rrna_table_filtered)

    # 2. Add informations from the classified fasta of MCHelper
    classified_module_fasta_file = mchelper_classified_fasta
    unclassified_module_fasta_file = mchelper_unclassified_fasta

    classified_module_records, classified_module_table = parse_fasta_lengths(classified_module_fasta_file)
    unclassified_module_records, unclassified_module_table = parse_fasta_lengths(unclassified_module_fasta_file)
    all_modules_table = pd.concat([classified_module_table, unclassified_module_table], ignore_index=True)
    te_table = merge_tables(te_table, all_modules_table)
    # fill missing values and ensure boolean dtype to avoid pandas future downcast warning
    te_table["replace_sequence"] = te_table["replace_sequence"].fillna(False).astype(bool)

    # 3. Add informations from the classification files of MCHelper
    classified_module_classif_file = mchelper_classif
    classif_table = parse_classif_file(classified_module_classif_file)
    te_table = merge_tables(te_table, classif_table, ["id", "comment"])

    # reorder columns to have "id" and "comment" first
    cols = te_table.columns.tolist()
    new_order = [cols[0], "comment"] + [c for c in cols if c not in (cols[0], "comment")]
    te_table = te_table[new_order]
    te_table["comment"] = te_table["comment"].fillna("")

    # Update rm_code using user-provided classification CSV when
    # MCHelper is present and rm_code == 'Unknown'. See user preference: only
    # update rows where MCHelper non-null AND rm_code == 'Unknown'.
    classification_df = parse_classification_csv(classification)
    te_table = update_rm_code_from_mchelper(te_table, classification_df)

    # 4. Look which sequences to filter out — apply filter and print removed rows
    te_table = filter_out_gene_rrna_hits(te_table)

    # 5. Save representative sequences for each family into a FASTA
    def save_representative_sequences(
        table: pd.DataFrame,
        classified_records: list,
        unclassified_records: list,
        fasta_path: Path,
        out_prefix: Path,
    ) -> None:
        """For each row in ``table``, choose a representative sequence and write to ``<out_prefix>.fa``.

        Rules:
        - If ``mchelper_id`` is present for the row, try to find the sequence in the
          MCHelper records (classified or unclassified) using the header base (before ``#``).
        - Otherwise, look for the sequence in the original FASTA prepared for MCHelper
          (``fasta_path``) by matching the header base (before ``#``) to the table ``id``.

        The output FASTA uses the table ``id`` as the sequence identifier.
        """

        # build lookup dictionaries keyed by the part before the first '#'
        mchelper_records = (classified_records or []) + (unclassified_records or [])
        mc_map: dict[str, SeqRecord] = {}
        for r in mchelper_records:
            key = str(r.id).split('#', 1)[0]
            mc_map[key] = r

        fasta_map: dict[str, SeqRecord] = {}
        for r in SeqIO.parse(fasta_path, 'fasta'):
            key = str(r.id).split('#', 1)[0]
            fasta_map[key] = r

        out_fasta = Path(out_prefix).with_suffix('.fa')
        out_fasta.parent.mkdir(parents=True, exist_ok=True)

        to_write: list[SeqRecord] = []
        missing_keys: list[str] = []

        for _, row in table.iterrows():
            # prefer mchelper_id if present
            m_id = row.get('mchelper_id')
            picked: SeqRecord | None = None

            if m_id is not None and str(m_id).strip() != '' and not pd.isna(m_id):
                key = str(m_id)
                picked = mc_map.get(key)
                if picked is None:
                    picked = fasta_map.get(key)
            else:
                key = str(row.get('id'))
                picked = fasta_map.get(key)

            # Replace the family part of the header (after '#') with the rm_code
            # that corresponds to this sequence's mchelper id in te_table.
            # Preserve any description text after the id.
            orig_id = str(picked.id)
            if '#' in orig_id:
                base, orig_family = orig_id.split('#', 1)
            else:
                base = orig_id
                orig_family = ''

            # find replacement rm_code using te_table mapping: prefer lookup by mchelper_id, fallback to id
            rm_code = None
            if base in te_table.get('mchelper_id', pd.Series(dtype='string')).astype(str).values:
                # select first matching
                sel = te_table.loc[te_table['mchelper_id'].astype(str) == base, 'rm_code']
                if not sel.empty:
                    rm_code = str(sel.iloc[0])
            if rm_code is None and base in te_table.get('id', pd.Series(dtype='string')).astype(str).values:
                sel = te_table.loc[te_table['id'].astype(str) == base, 'rm_code']
                if not sel.empty:
                    rm_code = str(sel.iloc[0])

            new_family = rm_code if rm_code not in (None, 'nan', '') else orig_family

            # build new id and preserve description by replacing the id occurrence
            new_id = f"{base}#{new_family}" if new_family != '' else base
            new_desc = picked.description.replace(orig_id, new_id, 1) if picked.description else new_id

            new_rec = SeqRecord(picked.seq, id=new_id, description=new_desc)
            # preserve annotations if present
            try:
                new_rec.annotations = picked.annotations
            except Exception:
                pass

            to_write.append(new_rec)

        if to_write:
            with out_fasta.open('w') as fh:
                SeqIO.write(to_write, fh, 'fasta')
            print(f"Wrote representative sequences to: {out_fasta}", file=sys.stderr)

        if missing_keys:
            print(f"Warning: {len(missing_keys)} families had no sequence found (examples: {missing_keys[:5]})", file=sys.stderr)

    # call sequence extraction
    try:
        save_representative_sequences(
            te_table,
            classified_module_records,
            unclassified_module_records,
            fasta_file,
            out_prefix or Path('.' ) / 'postprocess_mchelper',
        )
    except Exception as exc:
        print(f"Warning: failed to write representative sequences: {exc}", file=sys.stderr)

    # write final merged table if requested
    out_path = Path(out_prefix).with_suffix('.tsv')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    te_table.to_csv(out_path, sep='\t', index=False)
    print(f"Wrote merged table to: {out_path}", file=sys.stderr)

if __name__ == '__main__':
    main()
