// Include local modules
include { BLASTP                          } from "$projectDir/modules/local/module_genes/functional_annotation/blastp.nf"
include { INTERPROSCAN                    } from "$projectDir/modules/local/module_genes/functional_annotation/interproscan.nf"
include { EGGNOG_MAPPER                   } from "$projectDir/modules/local/module_genes/functional_annotation/eggnog_mapper.nf"
include { TRANSFER_FUNCTIONAL_ANNOTATIONS } from "$projectDir/modules/local/module_genes/functional_annotation/transfer_functional_annotations.nf"


workflow FUNCTIONAL_ANNOTATION {

    take:
        /*
         * Filtered structural GFF annotations.
         *
         * Expected format:
         *   [id, meta, annotation_method, gff]
         */
        gff

        /*
         * Canonical protein FASTA files extracted from the filtered GFF annotations.
         *
         * Expected format:
         *   [id, meta, annotation_method, fasta]
         */
        proteins

    main:

        /*
         * Prepare BLASTP input against Swiss-Prot/UniProt.
         *
         * Output format:
         *   [id, meta, annotation_method, fasta, db_name, evalue, max_target_seqs, num_descriptions]
         */
        proteins
            .map { id, meta, annotation_method, fasta ->
                tuple(
                    id,
                    meta,
                    annotation_method,
                    fasta,
                    "uniprot_sprot",
                    "1e-6",
                    "3",
                    "20"
                )
            }
            .set { blastp_input }

        BLASTP(blastp_input)


        /*
         * Run InterProScan on canonical protein sequences.
         *
         * Expected input:
         *   [id, meta, annotation_method, fasta]
         */
        INTERPROSCAN(proteins)


        /*
         * Run eggNOG-mapper on canonical protein sequences.
         *
         * Expected input:
         *   [id, meta, annotation_method, fasta]
         */
        EGGNOG_MAPPER(proteins)


        /*
         * Combine structural GFF annotations with functional annotation results.
         *
         * Join key:
         *   [id, meta, annotation_method]
         *
         * Input to TRANSFER_FUNCTIONAL_ANNOTATIONS:
         *   [id, meta, annotation_method, gff, blastp_tsv, interproscan_tsv, eggnog_annotations]
         */
        gff
            .join(BLASTP.out.tsv, by: [0, 1, 2])
            .join(INTERPROSCAN.out.tsv, by: [0, 1, 2])
            .join(EGGNOG_MAPPER.out.annotations, by: [0, 1, 2])
            .set { transfer_functional_annotations_input }

        TRANSFER_FUNCTIONAL_ANNOTATIONS(transfer_functional_annotations_input)


    emit:
        /*
         * GFF files enriched with functional annotations.
         *
         * Output format:
         *   [id, meta, annotation_method, gff]
         */
        gff = TRANSFER_FUNCTIONAL_ANNOTATIONS.out

        /*
         * Raw functional annotation outputs.
         */
        blastp       = BLASTP.out.tsv
        interproscan = INTERPROSCAN.out.tsv
        eggnog       = EGGNOG_MAPPER.out.annotations
}