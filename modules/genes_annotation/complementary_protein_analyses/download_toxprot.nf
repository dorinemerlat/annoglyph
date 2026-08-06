process DOWNLOAD_TOXPROT {

    tag "toxprot"
    label 'exogap_download'
    scratch false

    output:
    path("blastdb"), emit: blastdb
    path("toxprot.fasta"), emit: fasta

    script:
    """
    wget -O toxprot.fasta.gz \\
        "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28taxonomy_id%3A33208%29+AND+%28%28cc_tissue_specificity%3Avenom+OR+cc_scl_term%3ASL-0177%29+OR+%28keyword%3AKW-0800%29%29+AND+%28reviewed%3Atrue%29%29"

    gunzip -c toxprot.fasta.gz > toxprot.fasta

    mkdir -p blastdb

    makeblastdb \\
        -in toxprot.fasta \\
        -dbtype prot \\
        -parse_seqids \\
        -out blastdb/toxprot
    """

    stub:
    """
    touch toxprot.fasta
    mkdir -p blastdb
    touch blastdb/toxprot.phr
    touch blastdb/toxprot.pin
    touch blastdb/toxprot.psq
    """
}