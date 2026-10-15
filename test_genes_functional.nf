include { FUNCTIONAL_ANNOTATION } from './subworkflows/genes_annotation/functional_annotation.nf'
workflow {
    meta = [id:'smoke', annotation_method:'braker2']
    proteins = Channel.of([meta, file("${projectDir}/test/validation/genes_real.fasta", checkIfExists: true)])
    gff = Channel.of([meta, file("${projectDir}/test/test_full2/cache/genes_annotation/structural_annotation/agat_spkeeplongestisoform/braker2/saccharomyces-paradoxus-asm207905v1/saccharomyces-paradoxus-asm207905v1_braker2_canonical.gff", checkIfExists: true)])
    FUNCTIONAL_ANNOTATION(gff, proteins)
    FUNCTIONAL_ANNOTATION.out.gff.view { meta, gff ->
        assert gff.size() > 0
        "FUNCTIONAL_REAL_OK: ${gff.name} (${gff.size()} bytes)"
    }
}
