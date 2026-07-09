include { BRAKER2                    } from "$projectDir/modules/local/module_genes/structural_annotation/braker2.nf"
include { BRAKER3                    } from "$projectDir/modules/local/module_genes/structural_annotation/braker3.nf"
include { GALBA                      } from "$projectDir/modules/local/module_genes/structural_annotation/galba.nf"
include { REFORMAT_GFF3              } from "$projectDir/modules/local/module_genes/structural_annotation/reformat_gff3.nf"
include { AGAT_FILTER_GENE_BY_LENGTH } from "$projectDir/modules/local/module_genes/structural_annotation/agat_filter_gene_by_length.nf"
include { EXTRACT_CANONICAL_PROTEINS } from "$projectDir/modules/local/module_genes/structural_annotation/extract_canonical_proteins.nf"


workflow STRUCTURAL_ANNOTATION {

    take:
        masked_genomes
        bams

    main:

        /*
         * BRAKER3: genome + RNA-seq BAM + protein hints
         *
         * Input expected:
         *   masked_genomes: [id, meta, genome]
         *   bams:           [id, meta, aligner, bam]
         */
        masked_genomes
            .combine(bams, by: [0, 1])
            .map { id, meta, genome, aligner, bam ->
                [ id, meta, aligner, genome, bam, file(params.protein_set) ]
            }
            .set { braker3_input }

        BRAKER3(braker3_input)


        /*
         * BRAKER2 and GALBA: genome + protein hints
         */
        masked_genomes
            .map { id, meta, genome ->
                [ id, meta, genome, file(params.protein_set) ]
            }
            .set { braker2_and_galba_input }

        BRAKER2(braker2_and_galba_input)
        GALBA(braker2_and_galba_input)


        /*
         * Merge all structural predictions into one standard channel:
         *
         *   [id, meta, annotation_method, gff]
         */
        BRAKER3.out.gff
            .map { id, meta, rnaseq_aligner, gff ->
                [ id, meta, "braker-with-${rnaseq_aligner}", gff ]
            }
            .concat(
                BRAKER2.out.gff.map { id, meta, gff ->
                    [ id, meta, "braker2", gff ]
                }
            )
            .concat(
                GALBA.out.gff.map { id, meta, gff ->
                    [ id, meta, "galba", gff ]
                }
            )
            .filter { id, meta, annotation_method, gff ->
                gff.size() > 0
            }
            .set { raw_gff }


        /*
         * Standardize GFF3 format
         */
        REFORMAT_GFF3(raw_gff)


        /*
         * Remove very short genes / proteins.
         * Same threshold as in the current workflow: 90 aa.
         */
        REFORMAT_GFF3.out
            .map { id, meta, annotation_method, gff ->
                [ id, meta, annotation_method, gff, 90 ]
            }
            .set { gff_for_length_filtering }

        AGAT_FILTER_GENE_BY_LENGTH(gff_for_length_filtering)


        /*
         * Extract canonical proteins from filtered GFF.
         *
         * Input expected by EXTRACT_CANONICAL_PROTEINS:
         *   [id, meta, annotation_method, genome, gff]
         */
        AGAT_FILTER_GENE_BY_LENGTH.out
            .combine(masked_genomes, by: [0, 1])
            .map { id, meta, annotation_method, gff, genome ->
                [ id, meta, annotation_method, genome, gff ]
            }
            .set { gff_and_genome }

        EXTRACT_CANONICAL_PROTEINS(gff_and_genome)


    emit:
        raw_gff         = raw_gff
        reformatted_gff = REFORMAT_GFF3.out
        filtered_gff    = AGAT_FILTER_GENE_BY_LENGTH.out
        gff             = EXTRACT_CANONICAL_PROTEINS.out.gff
        proteins        = EXTRACT_CANONICAL_PROTEINS.out.fasta
}