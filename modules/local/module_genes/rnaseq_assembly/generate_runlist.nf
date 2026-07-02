process GENERATE_RUNLIST {
    tag "${id}"
    scratch false
    maxForks 1
    label "retry_with_backoff"

    input:
    tuple val(id), val(meta), val(specie_name)

    output:
    tuple val(id), val(meta), path("Runlist.txt")

    script:
    def sra_list = meta.SRA ? meta.SRA.replaceAll(';', ' ') : null
    def name    = specie_name.replaceAll(/\s+/, '_').replace('.','')
    def toks    = name.tokenize('_')
    def genus   = toks.size() > 0 ? toks[0] : name
    def species = toks.size() > 1 ? toks[1..-1].join('_') : ''
    
    if (sra_list) {
        """
        generate_runlist.py -g "${genus}" -s "${species}" -i "${sra_list}" -o Runlist.txt
        """
    } else {
        """
        generate_runlist.py -g "${genus}" -s "${species}" -o Runlist.txt
        """
    }

    stub:
    """
    touch Runlist.txt
    """
}