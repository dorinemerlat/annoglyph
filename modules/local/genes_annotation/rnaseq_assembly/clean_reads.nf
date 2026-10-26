process CLEAN_RNASEQ_READS {
    tag "${meta.genome_id}/${meta.run_id}/${meta.readset_signature}"
    cpus 4
    memory 4.GB
    maxForks 4
    scratch false
    stageInMode 'symlink'
    container "${workflow.containerEngine in ['singularity', 'apptainer'] ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d0/d013aad5427d824afe472e6607ea47685ff0181f1fb09e52a179e0ec39e43e88/data' : 'community.wave.seqera.io/library/fastp:1.3.6--4df8d6c11b471bde'}"
    publishDir "${params.outdir}/rnaseq/qc/${meta.genome_id}/${meta.run_id}/${meta.readset_signature}",
        mode: params.publish_dir_mode, pattern: '*.fastp.{json,html}', overwrite: true

    input:
    tuple val(meta), path(reads, stageAs: 'input*/*')

    output:
    tuple val(meta), path('*.clean.fastq.gz'), path('*.fastp.json'), emit: reads
    tuple val(meta), path('*.fastp.html'), emit: html
    path 'versions.yml', emit: versions

    script:
    def files = reads instanceof List ? reads : [reads]
    def mate = meta.single_end ? '' : "--in2 '${files[1]}' --out2 '${meta.id}_R2.clean.fastq.gz'"
    def clean = params.rnaseq_clean_reads
    def adapters = !clean || !params.rnaseq_trim_adapters ? '--disable_adapter_trimming' : (meta.single_end ? '' : '--detect_adapter_for_pe')
    def quality = clean ? "--qualified_quality_phred ${params.rnaseq_qualified_quality_phred} --unqualified_percent_limit ${params.rnaseq_unqualified_percent_limit} --n_base_limit ${params.rnaseq_max_n_bases} --length_required ${params.rnaseq_min_read_length}" : '--disable_quality_filtering --disable_length_filtering'
    def polyG = clean && params.rnaseq_trim_poly_g ? '--trim_poly_g' : '--disable_trim_poly_g'
    def outputs = meta.single_end ? "'${meta.id}_R1.clean.fastq.gz'" : "'${meta.id}_R1.clean.fastq.gz' '${meta.id}_R2.clean.fastq.gz'"
    """
    fastp --in1 '${files[0]}' --out1 '${meta.id}_R1.clean.fastq.gz' ${mate} \
        --thread ${task.cpus} ${adapters} ${quality} ${polyG} \
        --dont_eval_duplication \
        --json '${meta.id}.fastp.json' --html '${meta.id}.fastp.html'
    # fastp may create zero-byte files when no reads pass. Emit valid empty gzip
    # streams so reports/outputs remain readable even though alignment is skipped.
    for fq in ${outputs}; do
        if [ ! -s "\$fq" ]; then gzip -c < /dev/null > "\$fq"; fi
    done
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
    END_VERSIONS
    """

    stub:
    def files = reads instanceof List ? reads : [reads]
    def mate = meta.single_end ? '' : "gzip -c '${files[1]}' > '${meta.id}_R2.clean.fastq.gz'"
    def count = task.ext.stub_read_count == null ? 100L : task.ext.stub_read_count as long
    """
    gzip -c '${files[0]}' > '${meta.id}_R1.clean.fastq.gz'
    ${mate}
    echo '{"summary":{"after_filtering":{"total_reads":${count},"total_bases":${count * 100}}}}' > '${meta.id}.fastp.json'
    touch '${meta.id}.fastp.html'
    echo 'fastp: stub' > versions.yml
    """
}
