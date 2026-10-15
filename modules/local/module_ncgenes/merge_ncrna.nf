process MERGE_NCRNA {
    tag "${meta.id}"
    label 'annoglyph_python'
    input:
    tuple val(meta), path(trnascan), path(barrnap_euk), path(barrnap_mito), path(infernal)
    path rfam
    output:
    tuple val(meta), path("${meta.id}_ncrna.gff3"), emit: gff
    tuple val(meta), path("${meta.id}_ncrna_filtering.tsv"), emit: filtering
    script:
    """
    merge_ncrna.py --trnascan ${trnascan} \
        --barrnap ${barrnap_euk} ${barrnap_mito} \
        --infernal ${infernal} --families ${rfam}/family.txt \
        --prefix '${meta.id}' --output ${meta.id}_ncrna.gff3 \
        --audit ${meta.id}_ncrna_filtering.tsv
    """
    stub:
    """
    echo '##gff-version 3' > ${meta.id}_ncrna.gff3
    echo 'model\tseqid\tstart\tend\tstrand\tdecision' > ${meta.id}_ncrna_filtering.tsv
    """
}
