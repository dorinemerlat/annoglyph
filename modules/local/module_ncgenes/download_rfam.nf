process DOWNLOAD_RFAM {
    tag "${release}"
    label 'exogap_python'
    input:
    val release
    output:
    path 'rfam', emit: models
    script:
    """
    download_rfam.py --release '${release}' --output rfam
    """
    stub:
    """
    mkdir rfam
    touch rfam/Rfam.cm rfam/Rfam.clanin rfam/family.txt
    """
}

process RFAM_CMPRESS {
    tag "${params.rfam_release}"
    label 'infernal'
    input:
    path rfam
    output:
    path 'rfam_indexed', emit: models
    script:
    """
    cp -rL ${rfam} rfam_indexed
    cmpress rfam_indexed/Rfam.cm
    """
    stub:
    """
    mkdir rfam_indexed
    touch rfam_indexed/Rfam.cm rfam_indexed/Rfam.clanin rfam_indexed/family.txt
    touch rfam_indexed/Rfam.cm.i1f rfam_indexed/Rfam.cm.i1i rfam_indexed/Rfam.cm.i1m rfam_indexed/Rfam.cm.i1p
    """
}
