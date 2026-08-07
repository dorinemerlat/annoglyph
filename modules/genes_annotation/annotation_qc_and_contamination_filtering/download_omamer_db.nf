process DOWNLOAD_OMAMER_DB {
    output:
    path("LUCA.h5")

    script:
    """
    wget \
        --tries=3 \
        --timeout=60 \
        -O LUCA.h5 \
        https://omabrowser.org/All/LUCA.h5
    """

    stub:
    """
    touch LUCA.h5
    """
}
