process DOWNLOAD_NEWICK {
    tag "newick"
    label 'exogap_python'

    input:
    val(taxids)

    output:
    path "tree.nwk"

    script:
    """
    download_newick.py  --outfile tree.nwk $taxids
    """

    stub:
    """
    tree stub.nwk
    """
}
