nextflow.enable.dsl=2

process PLOT_REPEATS {
	// label 'exogap_tools'
	input:
	tuple val(meta), path(gff_files), path(stats_files), path(newick)

	output:
	path "TE_overview.pdf"
	path "heatmap_load_log.pdf"
	path "genome_n50.tsv"

	when:
	meta.id == 'all'

	script:
	"""
	mkdir -p plots_input
	mv *gff plots_input/
	mv *stats plots_input/
	plot_repeats.R
	"""
}

