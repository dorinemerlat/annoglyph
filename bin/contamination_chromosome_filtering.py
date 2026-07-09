#!/usr/bin/env python3

import gffutils
from Bio import SeqIO
import os
import sys
import argparse
import re


def clean_id(x):
    """
    Remove only version suffixes such as .1, .2, .12.
    Does NOT modify IDs like THETU.g1.t1.
    """
    if x is None:
        return None
    return re.sub(r'\.\d+$', '', x)


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description=(
            "This script uses the OMArk contaminant list to identify corresponding "
            "areas in a genome and filter these proteins out."
        )
    )

    parser.add_argument('-i', '--omark_folder', type=str, required=True)
    parser.add_argument('-g', '--gff', type=str, required=True)
    parser.add_argument('-f', '--fasta', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)

    parser.add_argument(
        '-t', '--threshold',
        type=float,
        default=0.5,
        help="Threshold between 0 and 1."
    )

    parser.add_argument(
        '-n', '--minimal_gene_nr',
        type=int,
        default=5,
        help="Minimum number of contaminant genes in a stretch."
    )

    parser.add_argument(
        '-x', '--only_extrem',
        action='store_true',
        help="Only consider full-contig or contig-extremity stretches."
    )

    return parser


def get_contaminants(omark_folder):
    omark_categ = None

    for file in os.listdir(omark_folder):
        if file.endswith('.ump'):
            omark_categ = os.path.join(omark_folder, file)
            break

    if not omark_categ:
        sys.exit('An .ump file was not found in your folder')

    conta = False
    contaminant_ids = []

    with open(omark_categ) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith('>'):
                if line[1:].split('_')[0] == 'Contamination':
                    conta = True
                else:
                    conta = False
            elif conta:
                contaminant_ids.append(clean_id(line))

    return contaminant_ids


def get_position_conta(contaminant_list, gff):
    gene_to_prot = {}
    contaminant_positions = []
    all_positions = []

    contaminant_set = set(clean_id(x) for x in contaminant_list)

    db = gffutils.create_db(
        gff,
        ':memory:',
        merge_strategy="create_unique",
        keep_order=True
    )

    for t in db.features_of_type('gene', order_by='start'):
        gene_id = t.id.replace('gene:', '')
        gene_is_cont = False

        ordered_child = list(db.children(t, featuretype='CDS', order_by='start'))

        for child in ordered_child:
            protein = None

            for att_type in ['protein_id', 'Name']:
                protein = child.attributes.get(att_type, [None])[0]
                if protein:
                    break

            if protein:
                protein_clean = clean_id(protein)

                if protein not in gene_to_prot.get(gene_id, []):
                    prot_in_gene = gene_to_prot.get(gene_id, [])
                    prot_in_gene.append(protein)
                    gene_to_prot[gene_id] = prot_in_gene

                if protein_clean in contaminant_set:
                    gene_is_cont = True

        if gene_is_cont:
            contaminant_positions.append((t.chrom, t.start, t.end, gene_id))

        all_positions.append((t.chrom, t.start, t.end, gene_id))

    return contaminant_positions, all_positions, gene_to_prot


def infer_contaminant_genome_stretches(
    cont_pos,
    all_pos,
    thresh=0.5,
    min_number_in_stretch=5,
    force_extremities=True
):
    sorted_possible = sorted(all_pos, key=lambda element: (element[0], element[1]))
    cont_pos_set = set(cont_pos)

    stretch_list = []

    cur_start = None
    cur_end = None
    pos_in_stretch = 0
    gene_in_stretch = 0
    prev_chrom = None
    in_stretch = False
    pos_from_last = 0
    from_start = 9999
    from_last = 9999
    without_pos = 0

    for elem in sorted_possible:
        chrom, start, end, gene_id = elem

        if chrom != prev_chrom:
            if prev_chrom is not None and gene_in_stretch > 0:
                if (pos_in_stretch >= min_number_in_stretch) or (
                    from_start > 0 and pos_from_last / from_start == 1
                ):
                    if from_last > 0 and pos_from_last / from_last >= thresh:
                        gene_in_stretch += without_pos
                        fraction_gene = pos_in_stretch / gene_in_stretch
                        stretch_list.append(
                            (
                                prev_chrom,
                                cur_start,
                                -1,
                                pos_in_stretch,
                                gene_in_stretch,
                                fraction_gene
                            )
                        )
                    else:
                        fraction_gene = pos_in_stretch / gene_in_stretch
                        if not force_extremities or cur_end == -1 or cur_start == 1:
                            stretch_list.append(
                                (
                                    prev_chrom,
                                    cur_start,
                                    cur_end,
                                    pos_in_stretch,
                                    gene_in_stretch,
                                    fraction_gene
                                )
                            )

            prev_chrom = chrom
            from_start = 0
            from_last = 0
            without_pos = 0
            pos_from_last = 0
            pos_in_stretch = 0
            gene_in_stretch = 0
            cur_start = start
            cur_end = end
            in_stretch = False

        from_start += 1
        from_last += 1

        if elem in cont_pos_set:
            pos_from_last += 1

            if pos_from_last / from_start >= thresh:
                if not in_stretch:
                    cur_start = 1
                    cur_end = end
                    in_stretch = True
                    gene_in_stretch = from_start
                else:
                    cur_end = end
                    gene_in_stretch += without_pos + 1

                pos_in_stretch += 1

            elif pos_from_last / from_last >= thresh:
                cur_end = end
                pos_in_stretch += 1
                gene_in_stretch += without_pos + 1

            else:
                if in_stretch and pos_in_stretch >= min_number_in_stretch:
                    if not force_extremities or cur_end == -1 or cur_start == 1:
                        fraction_gene = pos_in_stretch / gene_in_stretch
                        stretch_list.append(
                            (
                                chrom,
                                cur_start,
                                cur_end,
                                pos_in_stretch,
                                gene_in_stretch,
                                fraction_gene
                            )
                        )

                pos_in_stretch = 1
                cur_start = start
                cur_end = end
                pos_from_last = 1
                from_last = 1
                gene_in_stretch = 1

            without_pos = 0

        else:
            without_pos += 1

    if prev_chrom is not None and gene_in_stretch > 0:
        if pos_in_stretch >= min_number_in_stretch:
            fraction_gene = pos_in_stretch / gene_in_stretch
            if not force_extremities or cur_end == -1 or cur_start == 1:
                stretch_list.append(
                    (
                        prev_chrom,
                        cur_start,
                        cur_end,
                        pos_in_stretch,
                        gene_in_stretch,
                        fraction_gene
                    )
                )

    return stretch_list


def get_genes_in_cont_stretches(contaminant_stretches, all_pos, gene_to_prot):
    prot_list = []
    gene_list = []

    sorted_possible = sorted(all_pos, key=lambda element: (element[0], element[1]))

    for elem in contaminant_stretches:
        for elem2 in sorted_possible:
            if elem2[0] < elem[0]:
                continue
            elif elem2[0] > elem[0]:
                break
            else:
                if elem2[1] >= elem[1] and (elem2[2] <= elem[2] or elem[2] == -1):
                    gene_list.append(elem2)
                    prot_list += gene_to_prot.get(elem2[3], [])

    return prot_list, gene_list


def filter_proteins(fasta, newfa, cont_proteins):
    keep_record = []

    cont_set = set(clean_id(x) for x in cont_proteins)

    total = 0
    removed = 0

    for record in SeqIO.parse(fasta, "fasta"):
        total += 1
        rec_id = clean_id(record.id)

        if rec_id not in cont_set:
            keep_record.append(record)
        else:
            removed += 1

    SeqIO.write(keep_record, newfa, "fasta")

    print(f"Total proteins in input FASTA: {total}")
    print(f"Proteins removed: {removed}")
    print(f"Proteins kept: {len(keep_record)}")


def write_report(
    contaminant_stretches,
    gene_to_remove,
    protein_to_remove,
    gene_to_prot,
    outfile
):
    all_str = []

    all_str.append('-- Contaminant contig segments --')
    all_str.append(
        '\t'.join(
            [
                'Chromosome',
                'Start',
                'End',
                'Number of contaminants genes',
                'Number of genes',
                'Fraction contaminant genes'
            ]
        )
    )

    for cont in contaminant_stretches:
        all_str.append('\t'.join([str(x) for x in cont]))

    all_str.append('-- Contaminant genes --')
    all_str.append('\t'.join(['Gene', 'Chromosome', 'Start', 'End']))

    for gene in gene_to_remove:
        all_str.append('\t'.join([gene[3]] + [str(x) for x in gene[0:3]]))

    all_str.append('-- Contaminant proteins --')
    all_str.append('\t'.join(['Protein', 'Gene']))

    prot_to_gene = {}

    for gene_id, prots in gene_to_prot.items():
        for prot in prots:
            prot_to_gene[prot] = gene_id

    for prot in protein_to_remove:
        all_str.append('\t'.join([prot, prot_to_gene.get(prot, "NA")]))

    with open(outfile, 'w') as of:
        of.write('\n'.join(all_str))


if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()

    omark_folder = arg.omark_folder
    gff_file = arg.gff
    original_fasta = arg.fasta
    output = arg.output

    threshold = arg.threshold
    min_num = arg.minimal_gene_nr
    only_extremities = arg.only_extrem

    outfasta = output + ".fa"
    outtext = output + ".txt"

    contaminant = get_contaminants(omark_folder)

    print(f"Contaminant proteins reported by OMArk: {len(contaminant)}")

    cont_pos, all_pos, gene_to_prot = get_position_conta(contaminant, gff_file)

    print(f"Genes in GFF: {len(all_pos)}")
    print(f"Genes with contaminant proteins: {len(cont_pos)}")

    contaminant_stretches = infer_contaminant_genome_stretches(
        cont_pos,
        all_pos,
        threshold,
        min_num,
        only_extremities
    )

    print(f"Contaminant stretches detected: {len(contaminant_stretches)}")

    protein_to_remove, gene_to_remove = get_genes_in_cont_stretches(
        contaminant_stretches,
        all_pos,
        gene_to_prot
    )

    print(f"Genes to remove: {len(gene_to_remove)}")
    print(f"Proteins to remove: {len(protein_to_remove)}")

    filter_proteins(original_fasta, outfasta, protein_to_remove)

    write_report(
        contaminant_stretches,
        gene_to_remove,
        protein_to_remove,
        gene_to_prot,
        outfile=outtext
    )