include { PREPARE_PROTEIN_EVIDENCE } from "$projectDir/modules/local/genes_annotation/structural_annotation/prepare_protein_evidence.nf"
include { BRAKER2 } from "$projectDir/modules/local/genes_annotation/structural_annotation/braker2.nf"
include { BRAKER3 } from "$projectDir/modules/local/genes_annotation/structural_annotation/braker3.nf"
include { GALBA } from "$projectDir/modules/local/genes_annotation/structural_annotation/galba.nf"
include { REFORMAT_GFF3 } from "$projectDir/modules/local/genes_annotation/structural_annotation/reformat_gff3.nf"
include { AGAT_FILTER_GENE_BY_LENGTH } from "$projectDir/modules/local/genes_annotation/structural_annotation/agat_filter_gene_by_length.nf"
include { AGAT_SPKEEPLONGESTISOFORM } from "$projectDir/modules/nf-core/agat/spkeeplongestisoform/main.nf"
include { AGAT_SPEXTRACTSEQUENCES } from "$projectDir/modules/nf-core/agat/spextractsequences/main.nf"

workflow STRUCTURAL_ANNOTATION {
    take:
    masked_genomes
    bams
    taxonomy_table
    busco_dataset

    main:
    evidence_metadata = Channel.empty()
    if (params.protein_set) {
        protein_evidence = Channel.value(file(params.protein_set, checkIfExists: true))
    } else {
        PREPARE_PROTEIN_EVIDENCE(taxonomy_table, busco_dataset)
        protein_evidence = PREPARE_PROTEIN_EVIDENCE.out.proteins.first()
        evidence_metadata = PREPARE_PROTEIN_EVIDENCE.out.metadata.map { path -> ['protein_evidence_sources', path] }.mix(PREPARE_PROTEIN_EVIDENCE.out.release.map { path -> ['uniprot_release', path] })
    }

    // Aligners carry different metadata; associate their BAMs by sample identifier.
    masked_genomes.map { meta, genome -> [meta.id, meta, genome] }
        .combine(bams.map { meta, bam -> [meta.id, meta.rnaseq_aligner, bam] }, by: 0)
        .combine(protein_evidence)
        .map { id, meta, genome, aligner, bam, proteins ->
            [meta + [rnaseq_aligner: aligner], genome, bam, proteins]
        }
        .set { braker3_inputs }
    BRAKER3(braker3_inputs)

    masked_genomes.combine(protein_evidence).set { protein_only_inputs }
    BRAKER2(protein_only_inputs)
    GALBA(protein_only_inputs)

    BRAKER3.out.gff
        .map { meta, gff -> [meta + [annotation_method: "braker-with-${meta.rnaseq_aligner}"], gff] }
        .concat(BRAKER2.out.gff.map { meta, gff -> [meta + [annotation_method: 'braker2'], gff] })
        .concat(GALBA.out.gff.map { meta, gff -> [meta + [annotation_method: 'galba'], gff] })
        .filter { meta, gff -> gff.size() > 0 }
        .set { raw_gff }
    REFORMAT_GFF3(raw_gff)
    AGAT_FILTER_GENE_BY_LENGTH(REFORMAT_GFF3.out.map { meta, gff -> [meta, gff, params.gene_min_length] })

    AGAT_SPKEEPLONGESTISOFORM(AGAT_FILTER_GENE_BY_LENGTH.out, [])
    // Keep the genome and annotation in one tuple: several methods share a genome.
    AGAT_SPKEEPLONGESTISOFORM.out.gff.map { meta, gff -> [meta.id, meta, gff] }
        .combine(masked_genomes.map { meta, genome -> [meta.id, genome] }, by: 0)
        .map { id, meta, gff, genome -> [meta, gff, genome] }
        .set { sequence_inputs }
    AGAT_SPEXTRACTSEQUENCES(sequence_inputs, [])

    emit:
    protein_evidence = protein_evidence
    evidence_metadata = evidence_metadata
    versions = AGAT_SPKEEPLONGESTISOFORM.out.versions.mix(AGAT_SPEXTRACTSEQUENCES.out.versions)
    raw_gff = raw_gff
    reformatted_gff = REFORMAT_GFF3.out
    filtered_gff = AGAT_FILTER_GENE_BY_LENGTH.out
    gff = AGAT_SPKEEPLONGESTISOFORM.out.gff
    proteins = AGAT_SPEXTRACTSEQUENCES.out.fasta
}
