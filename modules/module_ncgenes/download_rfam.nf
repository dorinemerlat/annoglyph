process DOWNLOAD_RFAM {
    tag "DOWNLOAD_RFAM_${id}"

    output:
    tuple path("family.txt"),                                               emit: families
    tuple path("Rfam.clanin"), path("Rfam.cm"), path("family_cleaned.txt"), emit: models

    script:
    """
    wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz
    gunzip family.txt.gz
    cat family.txt |cut -f 1,19 |sed "s/;\s*\$//g" |sed "s/; /,/g" > family_cleaned.txt

    wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
    wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
    gunzip Rfam.cm.gz
    """

    stub:
    """
    touch family.txt family_cleaned.txt Rfam.clanin Rfam.cm
    """
}
