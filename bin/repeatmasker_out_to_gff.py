#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys
import pandas as pd
import re

def parse_args():
    p = argparse.ArgumentParser(
        description="Convert RepeatMasker .out to GFF with 9 columns"
    )
    p.add_argument('-r', '--repeatmasker', required=True, help='RepeatMasker .out file')
    p.add_argument('-m', '--mchelper', required=False, help='mchelper output (optional mapping file)')
    p.add_argument('-o', '--output', required=True, help='Output GFF file')
    p.add_argument('-c', '--classification', required=True, help='Classification file')
    return p.parse_args()


def load_mchelper_map(path: str):
    """Load MCHelper TSV as a dict key -> row-dict.
    Index by multiple identifiers so lookups succeed:
      - id (e.g. RHYIM_TE196)
      - mchelper_id or MCHelper_id (e.g. RHYIM_TE196_inc)
      - rm2_id when present
    Also expects annotation columns: rm_code, rm2_class, rm2_order, rm2_superfamily, rm2_family, rm2_subfamily, rm2_subfamily2.
    Missing values are handled later as NA.
    """
    if not path:
        return {}
    mp: dict[str, dict] = {}
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
            header = None
            for raw in fh:
                line = raw.rstrip('\n')
                if not line:
                    continue
                parts = line.split('\t')
                # first non-empty line that has many columns -> header
                if header is None:
                    header = parts
                    continue
                # pad to header length
                if len(parts) < len(header):
                    parts = parts + [''] * (len(header) - len(parts))
                row = {header[i]: parts[i] for i in range(len(header))}
                # add using multiple possible keys
                keys = set()
                for k in ('id', 'mchelper_id', 'MCHelper_id', 'rm2_id'):
                    v = row.get(k)
                    if v:
                        keys.add(v)
                for k in keys:
                    mp[k] = row
    except Exception:
        # best-effort: ignore parsing errors
        pass
    return mp


def parse_rm_line(tokens):
    """Parse tokens from a RepeatMasker .out line (split by whitespace).
    Returns dict with needed fields or None if malformed.
    Expected columns (simplified):
      0: score
      1: perc div.
      2: perc del.
      3: perc ins.
      4: sequence (seqid)
      5: begin
      6: end
      7: (left)
      8: strand (+ or C)
      9: repeat name
      10: class/family
    """
    try:
        score = tokens[0]
        seqid = tokens[4]
        start = int(tokens[5])
        end = int(tokens[6])
        strand_tok = tokens[8]
        repeat_name = tokens[9]
        classfamily = tokens[10] if len(tokens) > 10 else None
        # positions in repeat (begin, end) may include parentheses for C strand; strip them
        rep_begin = tokens[11] if len(tokens) > 11 else ''
        rep_end = tokens[12] if len(tokens) > 12 else ''
        try:
            rep_begin_i = int(rep_begin.strip('()')) if rep_begin else None
        except Exception:
            rep_begin_i = None
        try:
            rep_end_i = int(rep_end.strip('()')) if rep_end else None
        except Exception:
            rep_end_i = None

        strand = '+' if strand_tok == '+' else ('-' if strand_tok.upper() == 'C' else '.')
        return {
            'score': score,
            'seqid': seqid,
            'start': start,
            'end': end,
            'strand': strand,
            'repeat_name': repeat_name,
            'classfamily': classfamily,
            'rep_begin': rep_begin_i,
            'rep_end': rep_end_i,
        }
    except Exception:
        return None

def create_mcp(
    id : str = '',
    type_: str = '',
    comment: str = '',
    class_: str = '',
    order: str = '',
    superfamily: str = '',
    family: str = '',
    subfamily: str = '',
    subfamily2: str = ''
):
    mcp = {
        'id': id,
        'comment': comment,
        'rm2_id': '',
        'rm2_length': '',
        'rm_code': '',
        'rm2_type': type_,
        'rm2_class': class_,
        'rm2_order': order,
        'rm2_superfamily': superfamily,
        'rm2_family': family,
        'rm2_subfamily': subfamily,
        'rm2_subfamily2': subfamily2,
        'gene_target': '',
        'gene_evalue': '',
        'gene_score': '',
        'rrna_target': '',
        'rrna_evalue': '',
        'rrna_score': '',
        'mchelper_id': '',
        'MCHelper': '',
        'extended_length': '',
        'replace_sequence': '',
        'MCHelper_length': '',
        'MCHelper_strand': '',
        'MCHelper_confused': '',
        'MCHelper_class': '',
        'MCHelper_order': '',
        'MCHelper_Wcode': '',
        'MCHelper_sFamily': '',
        'MCHelper_CI': '',
        'MCHelper_coding': '',
        'MCHelper_struct': '',
        'MCHelper_other': '',
        'MCHelper_Reason': '',
        'statut': ''
    }

    return mcp

def normalize(s: str) -> str:
    if not s:
        return ''
    return s.lower().replace('-', '').replace('_', '').strip()


def extract_classification_dict(classif_row):
    """
    Build a standardized classification dictionary from a classif_df row
    """
    return {
        'type' : classif_row.get('Type', ''),
        'class': classif_row.get('Class', ''),
        'class': classif_row.get('Class', ''),
        'order': classif_row.get('Order', ''),
        'superfamily': classif_row.get('SuperFamily', ''),
        'family': classif_row.get('Family', ''),
        'subfamily': classif_row.get('SubFamily', ''),
        'subfamily2': classif_row.get('Subfamily2', '')
    }

# to delete later
def load_unknown_identified(path):
    identified = {}

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or "#" not in line:
                continue

            parts = line.split("#", 1)

            key = parts[0].lstrip(">").strip()
            value = parts[1].strip()

            identified[key] = value

    return identified


def check_problematic_repeat_name(repeat_name, mcp, classif_df):
    if repeat_name in mcp:
        if mcp[repeat_name]['rm2_class'] == 'ClassI' or mcp[repeat_name]['rm2_class'] == 'ClassII':
            mcp[repeat_name]['rm2_type'] = 'Transposable Element'
        elif mcp[repeat_name]['rm2_class'] == 'Satellite' or mcp[repeat_name]['rm2_class'] == 'Simple':
            mcp[repeat_name]['rm2_type'] = 'Tandem repeat'
        elif mcp[repeat_name]['rm2_class'] == 'RNA':
            mcp[repeat_name]['rm2_type'] = 'Pseudogene'
        elif mcp[repeat_name]['rm_code'].lower() == 'unknown':
            mcp[repeat_name]['rm2_type'] = 'Transposable Element'
        else:
            mcp[repeat_name]['rm2_type'] = 'Transposable Element'
        return mcp[repeat_name]

    elif re.search(r"\(([ACGT]+)\)", repeat_name, flags=re.IGNORECASE):
        mcp_line = create_mcp( id=repeat_name, comment="from RepBase", type_='Tandem repeat' )
        return mcp_line

    elif "-rich" in repeat_name.lower() or "_rich" in repeat_name.lower():
        mcp_line = create_mcp( id=repeat_name, comment="from RepBase", type_='Low complexity' )
        return mcp_line

    elif ("-" in repeat_name or "_" in repeat_name) and not '_TE' in repeat_name:
        first_part = re.split(r"[-_]", repeat_name, maxsplit=1)[0]
        first_part_norm = normalize(first_part)

        if first_part_norm:
            cols = ['Order', 'SuperFamily', 'Family', 'SubFamily', 'Subfamily2']
            cols = [c for c in cols if c in classif_df.columns]

            if cols:
                norm_df = classif_df[cols].apply(
                    lambda col: (
                        col.astype(str)
                           .str.lower()
                           .str.replace('-', '', regex=False)
                           .str.replace('_', '', regex=False)
                           .str.strip()
                    )
                )

                matches = norm_df.eq(first_part_norm)
                if matches.any(axis=None):
                    row_idx, col_name = matches.stack().loc[lambda x: x].index[0]
                    classif_row = classif_df.loc[row_idx]

                    classif_dict = extract_classification_dict(classif_row)
                    mcp_line = create_mcp(
                        id=repeat_name,
                        comment="from RepBase",
                        type_=classif_dict['type'],
                        class_=classif_dict['class'],
                        order=classif_dict['order'],
                        superfamily=classif_dict['superfamily'],
                        family=classif_dict['family'],
                        subfamily=classif_dict['subfamily'],
                        subfamily2=classif_dict['subfamily2']
                    )
                    return mcp_line
                else:
                    mcp_line = create_mcp( id=repeat_name, comment="from RepBase", type_='Unknown' )
                    return mcp_line

    elif not '_TE' in repeat_name:
        mcp_line = create_mcp( id=repeat_name, comment="from RepBase", type_='Unknown' )
        return mcp_line

    else:
        print(f"Warning: repeat name '{repeat_name}' not found in MCHelper map and does not match known patterns; classifying as unknown", file=sys.stderr)
        sys.exit(2)

def main():
    args = parse_args()
    in_path = Path(args.repeatmasker)
    if not in_path.exists():
        print(f"Error: input .out not found: {in_path}", file=sys.stderr)
        sys.exit(2)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # open classification file (tsv file) in a dataframe
    classification_path = Path(args.classification)
    if not classification_path.exists():
        print(f"Error: classification file not found: {classification_path}", file=sys.stderr)
        sys.exit(2)

    classif_df = pd.read_csv(classification_path, sep=",", header=0, dtype=str, keep_default_na=False)

    # load mchelper mapping if provided
    mcp = load_mchelper_map(args.mchelper) # mcp = dictionnary with infornations about each family annotated with MCHelper / RepeatModeler

    written = 0
    idx = 0

    with in_path.open('r', encoding='utf-8', errors='ignore') as fh, out_path.open('w', encoding='utf-8') as outfh:
        for line in fh:
            # skip header or empty lines
            if not line.strip() or line.startswith('SW') or line.startswith('score'):
                continue
            tokens = line.strip().split()
            # some .out files have leading spaces; ensure at least 11 tokens
            if len(tokens) < 11:
                continue
            parsed = parse_rm_line(tokens)
            if not parsed:
                continue

            idx += 1
            seqid = parsed['seqid']
            start = parsed['start']
            end = parsed['end']
            score = parsed['score']
            strand = parsed['strand']
            repeat_name = parsed['repeat_name']
            phase = '.' # phase is '.' for repeats
            total_length = end - start + 1

            # get the classification of mcp
            classification = check_problematic_repeat_name(repeat_name, mcp, classif_df)

            matching_repeat = repeat_name
            # Build ID as RE_<prefix>_<idx>
            prefix = seqid.split('_')[0] if '_' in seqid else seqid
            repeat_id = f"RE_{prefix}_{idx}"

            # if classification does not exist, assign 'NA' to all fields
            if not classification:
                classification = create_mcp( id=repeat_name, type_='Unknown' )
            rm2_type = classification.get('rm2_type')
            rm2_type = rm2_type.replace(' ', '_').lower().replace('tandem_repeat', 'simple_repeat')
            if rm2_type =='simple_repeat' or rm2_type == 'transposable_element' or rm2_type == 'unknown':
                feature = 'repetitive_elements'
            rm2_class = classification.get('rm2_class') or 'NA'
            if feature == 'repetitive_elements' and rm2_type == 'transposable_element' and rm2_class == 'NA' :
                rm2_type = 'unknown'
            rm2_order = classification.get('rm2_order') or 'NA'
            rm2_super = classification.get('rm2_superfamily') or 'NA'
            rm2_family = classification.get('rm2_family') or 'NA'
            rm2_subfam = classification.get('rm2_subfamily') or 'NA'
            rm2_subfam2 = classification.get('rm2_subfamily2') or 'NA'
            rm2_code = classification.get('rm_code') or 'NA'
            comment = classification.get('comment') or 'NA'
            MCHelper_coding = classification.get('MCHelper_coding') or 'NA'
            MCHelper_struct = classification.get('MCHelper_struct') or 'NA'
            MCHelper_other = classification.get('MCHelper_other') or 'NA'
            MCHelper_Reason = classification.get('MCHelper_Reason') or 'NA'
            status = classification.get('statut') or 'NA'

            attrs = [
                f"ID={repeat_id}",
                f"total_length={total_length}",
                f"matching_repeat={matching_repeat}",
                f"repeat_type={rm2_type}",
                f"repeat_class={rm2_class}",
                f"repeat_order={rm2_order}",
                f"repeat_superfamily={rm2_super}",
                f"repeat_family={rm2_family}",
                f"repeat_subfamily={rm2_subfam}",
                f"repeat_subfamily2={rm2_subfam2}",
                f"repeatmasker_code={rm2_code}",
                f"position_in_repeat_begin={parsed['rep_begin'] if parsed['rep_begin'] is not None else 'NA'}",
                f"position_in_repeat_end={parsed['rep_end'] if parsed['rep_end'] is not None else 'NA'}",
                f"comment={comment}",
                f"MCHelper_coding={MCHelper_coding}",
                f"MCHelper_struct={MCHelper_struct}",
                f"MCHelper_other={MCHelper_other}",
                f"MCHelper_Reason={MCHelper_Reason}",
                f"status={status}",
            ]
            # ensure trailing semicolon
            attr_str = ';'.join(attrs) + ';'

            # columns: seqid, source, type, start, end, score, strand, phase, attributes
            outfh.write('\t'.join([
                seqid,
                'RepeatMasker',
                feature,
                str(start),
                str(end),
                str(score),
                strand,
                phase,
                attr_str,
            ]) + '\n')
            written += 1

    print(f"✅ Wrote {written} records to {out_path}")


if __name__ == '__main__':
    main()
