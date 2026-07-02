/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SIGNALP         } from "$projectDir/modules/local/module_genes/complementary_protein_analyses/signalp.nf"
include { DOWNLOAD_TOXPROT  } from "$projectDir/modules/local/module_genes/complementary_protein_analyses/download_toxprot.nf"
include { TOXPROT_BLASTP  } from "$projectDir/modules/local/module_genes/complementary_protein_analyses/toxprot_blastp.nf"


workflow COMPLEMENTARY_PROTEIN_ANALYSES {

    take:
        /*
         * Best canonical protein FASTA files selected after annotation comparison.
         *
         * Expected format:
         *   [id, meta, fasta]
         */
        best_proteins

    main:

        /*
         * Predict secreted proteins.
         */
        SIGNALP(best_proteins)


        /*
        * Search proteins against a toxin/venom reference database.
        */
        DOWNLOAD_TOXPROT()

        best_proteins
            .combine(DOWNLOAD_TOXPROT.out.blastdb)
            .set { toxprot_blastp_input }

        TOXPROT_BLASTP(toxprot_blastp_input)


        /*
         * Optional future analyses:
         *
         * - Pfam / InterPro domain filtering
         * - venom-related keyword annotation
         * - cysteine-rich protein detection
         * - protein length / composition filters
         * - OrthoFinder / comparative genomics
         * - secretome + toxin intersection
         */


    emit:
        signalp  = SIGNALP.out.summary
        toxprot  = TOXPROT_BLASTP.out
}