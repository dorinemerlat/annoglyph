process REPEATMASKER_REPEATMASKER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.5--pl5321hdfd78af_0':
        'quay.io/biocontainers/repeatmasker:4.1.5--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(lib)

    output:
    tuple val(meta), path("${prefix}.masked"), emit: masked
    tuple val(meta), path("${prefix}.out")   , emit: out
    tuple val(meta), path("${prefix}.tbl")   , emit: tbl
    tuple val(meta), path("${prefix}.gff")   , emit: gff, optional: true
    tuple val(meta), path("${prefix}.align"), emit: align
    tuple val(meta), path("${prefix}.cat*"), emit: cat, optional: true
    tuple val(meta), path("${meta.id}_repeatmasker.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ''
    prefix      = task.ext.prefix   ?: "${meta.id}"
    def lib_arg = lib               ? "-lib $lib"   : ''

    def out_fasta    = fasta.name.endsWith('.gz') ? fasta.name[0..-4] : fasta.name
    def fasta_gz_cmd = fasta.name.endsWith('.gz') ? "gunzip -c ${fasta} > ${out_fasta}" : ""

    """
    ${fasta_gz_cmd}
    RepeatMasker \\
        $lib_arg \\
        -pa ${task.cpus} \\
        -dir ${prefix} \\
        ${args} \\
        ${out_fasta} 2>&1 | tee ${meta.id}_repeatmasker.log

    mv $prefix/${out_fasta}.align ${prefix}.align
    for cat_file in $prefix/${out_fasta}.cat*; do
        [ -e "\${cat_file}" ] || continue
        cat_suffix=\${cat_file##*.cat}
        mv "\${cat_file}" "${prefix}.cat\${cat_suffix}"
    done
    mv $prefix/${out_fasta}.masked  ${prefix}.masked
    mv $prefix/${out_fasta}.out     ${prefix}.out
    mv $prefix/${out_fasta}.tbl     ${prefix}.tbl
    mv $prefix/${out_fasta}.out.gff ${prefix}.gff       || echo "GFF is not produced"

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        repeatmasker: \$(RepeatMasker -v | sed 's/RepeatMasker version //1')
    END_VERSIONS
    """

    stub:
    prefix          = task.ext.prefix       ?: "${meta.id}"
    def args        = task.ext.args         ?: ''
    def touch_gff   = args.contains('-gff') ? "touch ${prefix}.gff" : ''

    """
    touch ${prefix}.masked
    touch ${prefix}.out
    touch ${prefix}.tbl
    touch ${prefix}.align
    touch ${prefix}.cat.gz
    touch ${meta.id}_repeatmasker.log
    $touch_gff

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        repeatmasker: \$(RepeatMasker -v | sed 's/RepeatMasker version //1')
    END_VERSIONS
    """
}
