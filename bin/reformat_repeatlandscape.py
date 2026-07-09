#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

LEVEL_COLS = ["Type", "Class", "Order", "SuperFamily", "Family", "SubFamily", "Subfamily2"]

ALLOWED_TYPES = {
    "unknown",
    "transposable_element",
    "simple_repeat",
    "pseudogene",
    "low_complexity",
}


def normalize_type(t: str) -> str:
    if t is None:
        return "unknown"

    t = str(t).strip().lower().replace(" ", "_").replace("-", "_")

    if t in ("transposableelement", "te"):
        return "transposable_element"
    if t in ("tandem_repeat", "tandemrepeat"):
        return "simple_repeat"
    if t in ("pseudo_gene",):
        return "pseudogene"
    if t in ("lowcomplexity",):
        return "low_complexity"

    return t if t in ALLOWED_TYPES else "unknown"


def build_code_map(classif_csv: Path):
    df = pd.read_csv(classif_csv, sep=",", dtype=str, keep_default_na=False)

    code_map = {}

    for _, row in df.iterrows():
        codes = [c.strip() for c in str(row["Code"]).split(";") if c.strip()]

        entry = {
            "Type": normalize_type(row.get("Type", "")),
            "Class": row.get("Class", "").strip() or "NA",
            "Order": row.get("Order", "").strip() or "NA",
            "SuperFamily": row.get("SuperFamily", "").strip() or "NA",
            "Family": row.get("Family", "").strip() or "NA",
            "SubFamily": row.get("SubFamily", "").strip() or "NA",
            "Subfamily2": row.get("Subfamily2", "").strip() or "NA",
        }

        for c in codes:
            code_map[c] = entry
            code_map[c.lower()] = entry

    return code_map


def entry_to_label(entry):
    return "|".join(entry[c] for c in LEVEL_COLS)


def parse_label(label):
    parts = str(label).split("|")
    if len(parts) == 7:
        return dict(zip(LEVEL_COLS, parts))
    return {
        "Type": "unknown",
        "Class": "NA",
        "Order": "NA",
        "SuperFamily": "NA",
        "Family": "NA",
        "SubFamily": "NA",
        "Subfamily2": "NA",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-c", "--classification", required=True)
    parser.add_argument("--id", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    in_path = Path(args.input)
    classif_path = Path(args.classification)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not in_path.exists():
        raise FileNotFoundError(in_path)

    code_map = build_code_map(classif_path)

    # ---------- read header ----------
    with in_path.open() as fh:
        header = fh.readline().rstrip("\n").split()

    div_col = header[0]
    codes = header[1:]

    # map codes to classification labels
    new_cols = []
    for c in codes:
        entry = code_map.get(c) or code_map.get(c.lower())
        if entry:
            new_cols.append(entry_to_label(entry))
        else:
            new_cols.append("unknown|NA|NA|NA|NA|NA|NA")

    # ---------- load dataframe ----------
    df = pd.read_csv(in_path, delim_whitespace=True)
    df.columns = [div_col] + new_cols

    # ---------- wide -> long ----------
    long_df = df.melt(
        id_vars=[div_col],
        var_name="repeat_label",
        value_name="percentage",
    )

    # drop zeros by default
    long_df = long_df[long_df["percentage"] != 0].copy()

    # IMPORTANT: réaligne les index
    long_df = long_df.reset_index(drop=True)

    levels = long_df["repeat_label"].apply(parse_label).apply(pd.Series)

    out_df = pd.concat(
        [
            pd.Series(args.id, index=long_df.index, name="genome_id"),
            levels,
            long_df[[div_col, "percentage"]].rename(columns={div_col: "Div"}),
        ],
        axis=1,
    )

    out_df = out_df.dropna(subset=["percentage"])

    out_df.columns = (
    out_df.columns
        .str.lower()
        .str.replace(" ", "_")
    )

    out_df.to_csv(out_path, sep="\t", index=False)

    print(f"Wrote {len(out_df)} rows → {out_path}")


if __name__ == "__main__":
    main()
