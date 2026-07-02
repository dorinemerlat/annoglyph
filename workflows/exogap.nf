/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPROCESSING } from "$projectDir/subworkflows/local/preprocessing.nf"
include { REPEATS_ANNOTATION } from "$projectDir/subworkflows/local/repeats_annotation.nf"
include { GENES_ANNOTATION } from "$projectDir/subworkflows/local/genes_annotation/genes_annotation.nf"
include { NCGENES_ANNOTATION } from "$projectDir/subworkflows/local/ncgenes_annotation.nf"
include { POSTPROCESSING } from "$projectDir/subworkflows/local/postprocessing.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EXOGAP {
	take:
		genomes

	main:
		log.info "Running workflow EXOGAP"

		log.info "Running module cleaning genomes"
		PREPROCESSING(genomes)

		if (params.module_repeats == true) {
			log.info "Running module repeat annotation"
			REPEATS_ANNOTATION(PREPROCESSING.out.genomes, PREPROCESSING.out.newick, PREPROCESSING.out.genome_stats)

		}

		if (params.module_genes == true) {
			log.info "Running module gene annotation"

			if (params.module_repeats == true) {
				REPEATS_ANNOTATION.out.masked_genomes.set { masked_genomes }
				REPEATS_ANNOTATION.out.repeat_gff.set { repeat_gff }
			} else {
				PREPROCESSING.out.genomes.set { masked_genomes }
				Channel.empty().set { repeat_gff }
			}

			GENES_ANNOTATION(
				masked_genomes,
				repeat_gff,
				PREPROCESSING.out.newick,
				PREPROCESSING.out.genome_stats
			)
		}

		// if (params.module_ncgenes == true) {
		//     print "Running module ncRNA annotation"

		//     if (params.module_repeats == true) {
		//         REPEATS_ANNOTATION.out.unmasked_genomes.set { genomes_for_ncgenes_annotation }
		//     } else {
		//         PREPROCESSING.out.genomes.set { genomes_for_ncgenes_annotation }
		//     }

		//     NCGENES_ANNOTATION( genomes_for_ncgenes_annotation )
		// }

		// print "Running module postprocessing"
		// POSTPROCESSING()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
	if (params.email || params.email_on_fail) {
		NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
	}
	NfcoreTemplate.summary(workflow, params, log)
	if (params.hook_url) {
		NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
	}
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
