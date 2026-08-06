include { PREPARE_PROTEIN_EVIDENCE   } from "$projectDir/modules/genes_annotation/structural_annotation/prepare_protein_evidence.nf"
include { BRAKER2                    } from "$projectDir/modules/genes_annotation/structural_annotation/braker2.nf"
include { BRAKER3                    } from "$projectDir/modules/genes_annotation/structural_annotation/braker3.nf"
include { GALBA                      } from "$projectDir/modules/genes_annotation/structural_annotation/galba.nf"
include { REFORMAT_GFF3              } from "$projectDir/modules/genes_annotation/structural_annotation/reformat_gff3.nf"
include { AGAT_FILTER_GENE_BY_LENGTH } from "$projectDir/modules/genes_annotation/structural_annotation/agat_filter_gene_by_length.nf"
include { EXTRACT_CANONICAL_PROTEINS } from "$projectDir/modules/genes_annotation/structural_annotation/extract_canonical_proteins.nf"


workflow STRUCTURAL_ANNOTATION {
    take:
        masked_genomes
        bams
        taxonomy_table
        busco_datasets_tsv

    main:
        // Use user-provided or automatically generated protein evidence.
        if (params.protein_set) {
            protein_evidence = Channel.value(file(params.protein_set, checkIfExists: true))
        } else {
            PREPARE_PROTEIN_EVIDENCE(taxonomy_table, busco_datasets_tsv)
            protein_evidence = PREPARE_PROTEIN_EVIDENCE.out.proteins
        }

        // Run BRAKER3 with RNA-seq and protein evidence.
        masked_genomes
            .combine(bams, by: [0, 1])
            .combine(protein_evidence)
            .map { id, meta, genome, aligner, bam, proteins -> [id, meta, aligner, genome, bam, proteins] }
            .set { braker3_inputs }

        BRAKER3(braker3_inputs)

        // Run BRAKER2 and GALBA with protein evidence.
        masked_genomes
            .combine(protein_evidence)
            .map { id, meta, genome, proteins -> [id, meta, genome, proteins] }
            .set { braker2_and_galba_inputs }

        BRAKER2(braker2_and_galba_inputs)
        GALBA(braker2_and_galba_inputs)


        // Combine structural annotations.
        BRAKER3.out.gff
            .map { id, meta, aligner, gff -> [id, meta, "braker-with-${aligner}", gff] }
            .concat(BRAKER2.out.gff.map { id, meta, gff -> [id, meta, 'braker2', gff] })
            .concat(GALBA.out.gff.map { id, meta, gff -> [id, meta, 'galba', gff] })
            .filter { id, meta, annotationMethod, gff -> gff.size() > 0 }
            .set { raw_gff }

        // Standardize GFF3 files.
        REFORMAT_GFF3(raw_gff)

        // Remove genes encoding proteins shorter than 90 amino acids.
        REFORMAT_GFF3.out
            .map { id, meta, annotationMethod, gff -> [id, meta, annotationMethod, gff, 90] }
            .set { gff_for_length_filtering }

        AGAT_FILTER_GENE_BY_LENGTH(gff_for_length_filtering)

        // Extract canonical proteins from filtered annotations.
        AGAT_FILTER_GENE_BY_LENGTH.out
            .combine(masked_genomes, by: [0, 1])
            .map { id, meta, annotationMethod, gff, genome -> [id, meta, annotationMethod, genome, gff] }
            .set { gff_and_genome }

        EXTRACT_CANONICAL_PROTEINS(gff_and_genome)

    emit:
        raw_gff         = raw_gff
        reformatted_gff = REFORMAT_GFF3.out
        filtered_gff    = AGAT_FILTER_GENE_BY_LENGTH.out
        gff             = EXTRACT_CANONICAL_PROTEINS.out.gff
        proteins        = EXTRACT_CANONICAL_PROTEINS.out.fasta
}