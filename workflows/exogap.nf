include { PREPROCESSING         } from "$projectDir/subworkflows/preprocessing.nf"
include { REPEATS_ANNOTATION    } from "$projectDir/subworkflows/repeats_annotation.nf"
include { GENES_ANNOTATION      } from "$projectDir/subworkflows/genes_annotation/genes_annotation.nf"
// include { NCGENES_ANNOTATION    } from "$projectDir/subworkflows/ncgenes_annotation.nf"
// include { POSTPROCESSING        } from "$projectDir/subworkflows/postprocessing.nf"


workflow EXOGAP {
    take:
        genomes

    main:

        PREPROCESSING(genomes)

        if (params.module_repeats == true) {

            REPEATS_ANNOTATION(
                PREPROCESSING.out.genomes, 
                PREPROCESSING.out.newick, 
                PREPROCESSING.out.genome_stats, 
                PREPROCESSING.out.busco_datasets
            )
        }

        if (params.module_genes == true) {

            if (params.module_repeats == true) {
                masked_genomes = REPEATS_ANNOTATION.out.masked_genomes
                repeat_gff = REPEATS_ANNOTATION.out.repeat_gff
            } else {
                masked_genomes = PREPROCESSING.out.genomes
                Channel.empty().set { repeat_gff }
            }
            
            GENES_ANNOTATION(
                PREPROCESSING.out.genomes,
                masked_genomes, 
                repeat_gff, 
                PREPROCESSING.out.newick,
                PREPROCESSING.out.genome_stats, 
                PREPROCESSING.out.taxonomy_table,
                PREPROCESSING.out.busco_dataset_table
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

