process DOWNLOAD_GO_DB {
    tag "go"

    output:
    path("go-basic.obo")
    path("goslim_generic.obo")

    script:
    """
    wget -O go-basic.obo "http://current.geneontology.org/ontology/go-basic.obo"
    wget -O goslim_generic.obo "http://current.geneontology.org/ontology/subsets/goslim_generic.obo"
    """

    stub:
    """
    touch go-basic.obo goslim_generic.obo
    """
}