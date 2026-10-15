#!/usr/bin/env python3
"""Combine selected annotations and export consistent GFF3/FASTA/GenBank inputs."""
import argparse
import csv
import hashlib
import json
import re
import shutil
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from submission_support import load_metadata, modifiers
from exogap_gff import Feature, namespace, read_gff, sort_key


class Genome:
    def __init__(self, path):
        self.index = SeqIO.index(path, 'fasta')
        self.lengths = {seqid: len(self[seqid]) for seqid in self.index}

    @lru_cache(maxsize=2)
    def __getitem__(self, seqid):
        return self.index[seqid]

    def __contains__(self, seqid):
        return seqid in self.index

    def __iter__(self):
        return iter(self.index)

    def close(self):
        self.__getitem__.cache_clear()
        self.index.close()


def write_fasta(handle, identifier, sequence, description=''):
    handle.write('>' + identifier + (' ' + description if description else '') + '\n')
    for offset in range(0, len(sequence), 60):
        handle.write(str(sequence[offset:offset + 60]) + '\n')


def descendants(identifier, children):
    pending, seen, result = list(children.get(identifier, [])), set(), []
    while pending:
        feature = pending.pop()
        key = (feature.attrs['ID'][0], feature.start, feature.end)
        if key in seen:
            continue
        seen.add(key)
        result.append(feature)
        pending.extend(children.get(feature.attrs['ID'][0], []))
    return result


def extract(parts, genome, strand):
    ordered = sorted(parts, key=lambda f: f.start, reverse=(strand == '-'))
    sequence = Seq('')
    for part in ordered:
        fragment = genome[part.seqid].seq[part.start - 1:part.end]
        sequence += fragment.reverse_complement() if strand == '-' else fragment
    return str(sequence).upper(), ordered


def qualifiers(handle, values):
    for key, value in values:
        if value:
            handle.write(f'\t\t\t{key}\t{value}\n')


def tbl_feature(handle, parts, kind, phase=0, partial_start=False, partial_end=False):
    for index, part in enumerate(parts):
        first, last = (part.end, part.start) if part.strand == '-' else (part.start, part.end)
        # GenBank partial markers follow biological 5' and 3' ends.
        left = ('>' if part.strand == '-' else '<') if partial_start and index == 0 else ''
        right = ('<' if part.strand == '-' else '>') if partial_end and index == len(parts) - 1 else ''
        handle.write(f'{left}{first}\t{right}{last}\t{kind if index == 0 else ""}\n')
    if kind == 'CDS':
        qualifiers(handle, [('codon_start', str(phase + 1))])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genome', required=True)
    parser.add_argument('--coding')
    parser.add_argument('--repeats')
    parser.add_argument('--ncrna')
    parser.add_argument('--sample', required=True)
    parser.add_argument('--organism', required=True)
    parser.add_argument('--taxid', required=True)
    parser.add_argument('--method', default='none')
    parser.add_argument('--locus-prefix', required=True)
    parser.add_argument('--translation-table', type=int, default=1)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--submission-metadata', type=Path)
    parser.add_argument('--assets', type=Path)
    args = parser.parse_args()
    submission, sequence_metadata = load_metadata(args.submission_metadata, args.sample, args.assets)
    args.locus_prefix = submission.get('locus_tag_prefix', args.locus_prefix)
    default_code = int(submission.get('translation_table', args.translation_table))
    if default_code not in CodonTable.unambiguous_dna_by_id:
        parser.error('Unknown genetic code')
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', args.locus_prefix):
        parser.error('Locus prefix must start with a letter and contain letters, digits or underscores')
    args.output.mkdir()
    out = args.output
    genome = Genome(args.genome)
    warnings, mapping_rows, features = [], [], []
    unknown_sequences = set(sequence_metadata) - set(genome.index)
    if unknown_sequences:
        raise ValueError(f'Sequence metadata do not match genome identifiers: {sorted(unknown_sequences)}')
    if default_code != args.translation_table:
        warnings.append([args.sample, 'submission_code_differs_from_prediction_code'])
    for source, path in [('coding', args.coding), ('repeats', args.repeats), ('ncrna', args.ncrna)]:
        if not path:
            continue
        records = read_gff(path)
        mapping = namespace(records, f'{args.sample}.{source}')
        mapping_rows.extend((source, old, new) for old, new in mapping.items())
        if source == 'repeats':
            for feature in records:
                feature.attrs.setdefault('original_type', [feature.kind])
                feature.kind = 'repeat_region'
        features.extend(records)
    # Some ncRNA predictors emit RNA rows without an enclosing gene.
    for feature in list(features):
        if (feature.kind.endswith('RNA') and feature.kind != 'mRNA'
                and 'Parent' not in feature.attrs):
            identifier = feature.attrs['ID'][0] + '.gene'
            features.append(Feature(feature.seqid, feature.source, 'gene', feature.start,
                                    feature.end, '.', feature.strand, '.', {'ID': [identifier]}))
            feature.attrs['Parent'] = [identifier]
    children = defaultdict(list)
    by_id = {}
    for feature in features:
        if feature.seqid not in genome or feature.end > genome.lengths[feature.seqid]:
            raise ValueError(f'Feature outside genome: {feature.line()}')
        by_id.setdefault(feature.attrs['ID'][0], feature)
        for parent in feature.attrs.get('Parent', []):
            children[parent].append(feature)
    for parent in children:
        if parent not in by_id:
            raise ValueError(f'Unresolved Parent: {parent}')
    for index, gene in enumerate(sorted((f for f in features if f.kind == 'gene'), key=sort_key), 1):
        gene.attrs['locus_tag'] = [f'{args.locus_prefix}_{index:06d}']

    def locus_for(feature):
        current, visited = feature, set()
        while current.kind != 'gene':
            parents = current.attrs.get('Parent', [])
            if not parents or parents[0] in visited:
                return ''
            visited.add(parents[0])
            current = by_id[parents[0]]
        return current.attrs['locus_tag'][0]

    cds_groups = []
    rna_groups = []
    with (out / 'cds.fna').open('w') as cds_handle, (out / 'proteins.faa').open('w') as prot_handle, \
            (out / 'ncrna.fna').open('w') as rna_handle:
        for transcript in sorted(features, key=sort_key):
            identifier = transcript.attrs['ID'][0]
            if transcript.kind == 'mRNA':
                parts = [f for f in descendants(identifier, children) if f.kind == 'CDS']
                if not parts:
                    warnings.append([identifier, 'mRNA_without_CDS'])
                    continue
                sequence, ordered = extract(parts, genome, transcript.strand)
                phase = int(ordered[0].phase) if ordered[0].phase in ('0', '1', '2') else 0
                table_values = transcript.attrs.get('transl_table', []) or ordered[0].attrs.get('transl_table', [])
                sequence_code = sequence_metadata.get(transcript.seqid, {}).get('translation_table')
                table = int(sequence_code) if sequence_code is not None else (int(table_values[0]) if table_values else default_code)
                if sequence_code is not None and table_values and int(table_values[0]) != table:
                    warnings.append([identifier, 'feature_code_differs_from_sequence_metadata'])
                if table not in CodonTable.unambiguous_dna_by_id:
                    raise ValueError(f'Unknown genetic code {table}')
                if sequence_metadata.get(transcript.seqid, {}).get('organelle'):
                    warnings.append([identifier, 'organelle_coding_prediction_requires_specific_review'])
                coding_sequence = sequence[phase:]
                usable = coding_sequence[:len(coding_sequence) // 3 * 3]
                translated = str(Seq(usable).translate(table=table))
                protein = translated[:-1] if translated.endswith('*') else translated
                product = transcript.attrs.get('product', ['hypothetical protein'])[0]
                product = re.sub(r'\s+', ' ', product)
                locus = locus_for(transcript)
                protein_id = identifier + '.protein'
                transcript.attrs['transcript_id'] = [f'gnl|EXOGAP|{identifier}']
                partial_start = phase > 0 or usable[:3] not in CodonTable.unambiguous_dna_by_id[table].start_codons
                if not partial_start and protein:
                    protein = 'M' + protein[1:]
                partial_end = len(coding_sequence) % 3 != 0 or not translated.endswith('*')
                cumulative = 0
                for i, part in enumerate(ordered):
                    expected = phase if i == 0 else (3 - cumulative % 3) % 3
                    if part.phase not in ('0', '1', '2') or int(part.phase) != expected:
                        warnings.append([identifier, f'inconsistent_CDS_phase:{part.start}-{part.end}'])
                    cumulative += part.end - part.start + 1 - (phase if i == 0 else 0)
                    part.attrs['protein_id'] = [f'gnl|EXOGAP|{protein_id}']
                    part.attrs['product'] = [product]
                    part.attrs['transl_table'] = [str(table)]
                    if locus:
                        part.attrs['locus_tag'] = [locus]
                if '*' in protein:
                    warnings.append([identifier, 'internal_stop_codon'])
                if partial_start or partial_end:
                    warnings.append([identifier, 'partial_CDS_requires_submission_review'])
                write_fasta(cds_handle, identifier, sequence, f'locus_tag={locus}')
                write_fasta(prot_handle, protein_id, protein,
                            f'{product} OS={args.organism} OX={args.taxid} GN={locus}')
                cds_groups.append((transcript, ordered, phase, product, table, partial_start, partial_end))
            elif transcript.kind.endswith('RNA'):
                related = descendants(identifier, children)
                parts = [f for f in related if f.kind == 'exon']
                if not parts:
                    # tRNAscan-SE may describe introns instead of explicit exons.
                    introns = sorted((f for f in related if f.kind == 'intron'), key=lambda f: f.start)
                    start = transcript.start
                    parts = []
                    for intron in introns:
                        if intron.start > start:
                            parts.append(Feature(transcript.seqid, transcript.source, 'exon', start, intron.start - 1, '.', transcript.strand, '.', {}))
                        start = intron.end + 1
                    if start <= transcript.end:
                        parts.append(Feature(transcript.seqid, transcript.source, 'exon', start, transcript.end, '.', transcript.strand, '.', {}))
                sequence, ordered = extract(parts, genome, transcript.strand)
                write_fasta(rna_handle, identifier, sequence,
                            f'type={transcript.kind} locus_tag={locus_for(transcript)}')
                rna_groups.append((transcript, ordered))
    if not cds_groups:
        warnings.append([args.sample, 'no_selected_protein_coding_annotation'])
    with (out / 'annotation.gff3').open('w') as handle:
        handle.write('##gff-version 3\n')
        for seqid in genome:
            handle.write(f'##sequence-region {seqid} 1 {genome.lengths[seqid]}\n')
        for feature in sorted(features, key=sort_key):
            handle.write(feature.line() + '\n')
    shutil.copyfile(args.genome, out / 'genome.masked.fna')
    # NCBI needs uppercase assembly sequence, independent of the analysis softmask.
    with (out / 'genome.fsa').open('w') as handle:
        for seqid in genome:
            write_fasta(handle, seqid, str(genome[seqid].seq).upper(),
                        modifiers(submission, sequence_metadata.get(seqid, {}), args.organism, default_code))
    with (out / 'annotation.tbl').open('w') as handle:
        genes_by_seq, cds_by_seq, rna_by_seq, repeats_by_seq = (defaultdict(list) for _ in range(4))
        for feature in features:
            if feature.kind == 'gene':
                genes_by_seq[feature.seqid].append(feature)
            elif feature.kind == 'repeat_region':
                repeats_by_seq[feature.seqid].append(feature)
        for group in cds_groups:
            cds_by_seq[group[0].seqid].append(group)
        for group in rna_groups:
            rna_by_seq[group[0].seqid].append(group)
        for seqid in genome:
            handle.write(f'>Feature {seqid}\n')
            for gene in sorted(genes_by_seq[seqid], key=sort_key):
                tbl_feature(handle, [gene], 'gene')
                qualifiers(handle, [('locus_tag', gene.attrs['locus_tag'][0])])
            for transcript, parts, phase, product, table, partial_start, partial_end in cds_by_seq[seqid]:
                identifier = transcript.attrs['ID'][0]
                exons = [f for f in descendants(identifier, children) if f.kind == 'exon']
                if exons:
                    ordered_exons = sorted(exons, key=lambda f: f.start, reverse=transcript.strand == '-')
                    tbl_feature(handle, ordered_exons, 'mRNA', partial_start=partial_start, partial_end=partial_end)
                    qualifiers(handle, [('product', product), ('transcript_id', f'gnl|EXOGAP|{identifier}')])
                else:
                    warnings.append([identifier, 'missing_exons_mRNA_omitted_from_feature_table'])
                tbl_feature(handle, parts, 'CDS', phase, partial_start, partial_end)
                qualifiers(handle, [('product', product), ('protein_id', f'gnl|EXOGAP|{identifier}.protein'),
                                    ('transcript_id', f'gnl|EXOGAP|{identifier}'), ('transl_table', str(table))])
            for transcript, parts in rna_by_seq[seqid]:
                kind = transcript.kind if transcript.kind in ('tRNA', 'rRNA') else 'ncRNA'
                tbl_feature(handle, parts, kind)
                product = transcript.attrs.get('product', transcript.attrs.get('Name', [kind]))[0]
                if kind == 'tRNA':
                    amino_acid = transcript.attrs.get('isotype', transcript.attrs.get('amino_acid', ['']))[0]
                    if not amino_acid:
                        match = re.search(r'(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Sec)', product)
                        amino_acid = match.group(1) if match else ''
                    if amino_acid:
                        product = 'tRNA-' + amino_acid
                    else:
                        warnings.append([transcript.attrs['ID'][0], 'tRNA_product_requires_review'])
                if kind == 'rRNA' and product.endswith('_rRNA'):
                    product = product[:-5] + ' ribosomal RNA'
                qualifiers(handle, [('product', product), ('locus_tag', locus_for(transcript))])
                if kind == 'ncRNA':
                    qualifiers(handle, [('ncRNA_class', 'other'), ('note', transcript.attrs.get('rfam_type', [transcript.kind])[0])])
            for repeat in repeats_by_seq[seqid]:
                tbl_feature(handle, [repeat], 'repeat_region')
                qualifiers(handle, [('note', repeat.attrs.get('repeat_family', ['RepeatMasker repeat'])[0])])
    with (out / 'id_mapping.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['source', 'original_id', 'export_id'])
        writer.writerows(mapping_rows)
    with (out / 'validation.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['feature', 'issue'])
        writer.writerows(warnings)
    metadata = {'sample': args.sample, 'organism': args.organism, 'taxid': args.taxid,
                'selected_method': args.method, 'locus_tag_prefix': args.locus_prefix,
                'coding_transcripts': len(cds_groups), 'noncoding_transcripts': len(rna_groups),
                'issues': len(warnings), 'submission_validated': False,
                'submission': submission, 'sequence_metadata': sequence_metadata, 'default_genetic_code': default_code,
                'ncRNA_annotation_enabled': bool(args.ncrna), 'repeat_annotation_enabled': bool(args.repeats),
                'required_before_submission': ['Registered BioProject/BioSample and locus_tag prefix',
                    'Submitter/source metadata and table2asn submission template',
                    'Resolve validation.tsv issues and run the destination validator']}
    (out / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    (out / 'README.txt').write_text(
        'One selected coding annotation plus repeat and noncoding annotations.\n'
        'ncrna.fna contains all ncRNA types, identified in each FASTA header.\n'
        'annotation.tbl and genome.fsa are inputs for NCBI table2asn, not a validated submission.\n'
        'The CDS FASTA preserves all annotated CDS bases; codon_start is represented in the feature table.\n'
        'Resolve validation.tsv and provide registered metadata/locus tags before submission.\n'
        'Protein FASTA identifiers are local IDs, not assigned UniProt accessions.\n')
    genome.close()


if __name__ == '__main__':
    main()
