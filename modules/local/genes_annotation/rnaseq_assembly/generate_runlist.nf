process GENERATE_RUNLIST {
    tag "${meta.id}"
    label 'exogap_python'

    input:
    val(meta)

    output:
    tuple val(meta), path("Runlist.txt"), emit: varus
    tuple val(meta), path("Runlist_filtered.txt"), emit: star
    tuple val(meta), path("run_metadata.tsv"), emit: metadata

    script:
    def scientificName = meta.lineage?.species?.name ?: meta.name
    def speciesParts = scientificName.tokenize(' ')
    def genus = speciesParts[0]
    def species = speciesParts.size() > 1 ? speciesParts[1] : ''
    def sraIds = meta.SRA ?: ''
    def runLimit = params.sra_run_limit == null || params.sra_run_limit instanceof Boolean ? 0 : params.sra_run_limit as int
    def maxBases = params.sra_run_max_bases == null ? 0 : params.sra_run_max_bases as long

    """
    generate_runlist.py \
        --genus '${genus}' \
        --species '${species}' \
        --sra-ids '${sraIds}' \
        --limit ${runLimit} \
        --max-bases ${maxBases} \
        --runlist-output Runlist.txt \
        --filtered-runlist-output Runlist_filtered.txt \
        --metadata-output run_metadata.tsv
    """

    stub:
    """
    command -v generate_runlist.py >/dev/null

    cat > Runlist.txt <<'EOF'
@Run_acc	total_spots	total_bases	avg_len	bool:paired	color_space
SRR000001	1000000	300000000	300.0	1	0
SRR000002	800000	120000000	150.0	0	0
EOF

    cat > Runlist_filtered.txt <<'EOF'
@Run_acc	total_spots	total_bases	avg_len	bool:paired	color_space
SRR000001	1000000	300000000	300.0	1	0
EOF

    cat > run_metadata.tsv <<'EOF'
Run_acc	total_spots	total_bases	avg_len	bool:paired	color_space	experiment_accession	study_accession	study_title	sample_accession	sample_name	experiment_title	library_name	library_strategy	library_source	library_selection	library_layout	instrument_model	platform
SRR000001	1000000	300000000	300.0	1	0	SRX000001	SRP000001	Test study	SRS000001	Test sample	Test experiment	Test library	RNA-Seq	TRANSCRIPTOMIC	cDNA	PAIRED	Illumina NovaSeq 6000	ILLUMINA
EOF
    """
}