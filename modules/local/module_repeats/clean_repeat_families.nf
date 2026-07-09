process CLEAN_REPEAT_FAMILIES {
    tag "${id}"
    label 'exogap_python'

    input:
    tuple val(id), val(meta), path(families)
    each path(classification)

    output:
    tuple val(id), val(meta), path("${families.baseName}.cleaned_for_mchelper.fa"), emit: for_mchelper
    tuple val(id), val(meta), path("${families.baseName}.tsv"), emit: table

    script:
    """
    clean_repeat_families.py -f ${families} -c ${classification} -m ${meta.mnemonic.exogap}
    """ 

    stub:
    """
    touch ${families.baseName}.cleaned.fa ${families.baseName}.cleaned_for_mchelper.fa
    """
}
