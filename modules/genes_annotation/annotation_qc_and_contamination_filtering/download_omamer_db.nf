process DOWNLOAD_OMAMER_DB {
    tag "luca"

    output:
    path("LUCA.h5")

    script:
    """
    wget "https://omabrowser.org/All/LUCA.h5"
    """

    stub:
    """
    touch LUCA.h5
    """
}
