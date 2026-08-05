nextflow.enable.dsl=2

process PLOT_REPEATS {
	tag "plot_repeats"
	// label 'exogap_tools'
	input:
	tuple val(id), path(gff_files), path(stats_files), path(newick)

	output:
	path "TE_overview.pdf"
	path "heatmap_load_log.pdf"
	path "genome_n50.tsv"

	when:
	id == 'all'

	script:
	"""
	mkdir -p plots_input
	mv *gff plots_input/
	mv *stats plots_input/
	plot_repeats.R
	"""
}

