include { BLASTP                          } from "$projectDir/modules/genes_annotation/functional_annotation/blastp.nf"
include { INTERPROSCAN                    } from "$projectDir/modules/genes_annotation/functional_annotation/interproscan.nf"
include { EGGNOG_MAPPER                   } from "$projectDir/modules/genes_annotation/functional_annotation/eggnog_mapper.nf"
include { TRANSFER_FUNCTIONAL_ANNOTATIONS } from "$projectDir/modules/genes_annotation/functional_annotation/transfer_functional_annotations.nf"
include { DOWNLOAD_TOXPROT                } from "$projectDir/modules/genes_annotation/functional_annotation/download_toxprot.nf"
include { TOXPROT_BLASTP                  } from "$projectDir/modules/genes_annotation/functional_annotation/toxprot_blastp.nf"
include { SIGNALP                         } from "$projectDir/modules/genes_annotation/functional_annotation/signalp.nf"


workflow FUNCTIONAL_ANNOTATION {

    take:
        gff
        proteins

    main:
        // Search canonical proteins against Swiss-Prot.
        proteins
            .map { id, meta, annotation_method, fasta ->
                [id, meta, annotation_method, fasta, 'uniprot_sprot', 1e-6, 3, 20]
            }
            .set { blastp_inputs }

        BLASTP(blastp_inputs)

        // Identify protein domains and functional signatures.
        INTERPROSCAN(proteins)

        // Assign orthology-based functional annotations.
        EGGNOG_MAPPER(proteins)

        // Transfer functional annotations to the GFF3.
        gff
            .join(BLASTP.out.tsv, by: [0, 1, 2])
            .join(INTERPROSCAN.out.tsv, by: [0, 1, 2])
            .join(EGGNOG_MAPPER.out.annotations, by: [0, 1, 2])
            .set { functional_annotation_inputs }

        TRANSFER_FUNCTIONAL_ANNOTATIONS(functional_annotation_inputs)

        // Detect signal peptides.
        if (params.signalp) {
            SIGNALP(proteins)
        }

        // Search proteins against ToxProt.
        if (params.toxprot_blast) {
            DOWNLOAD_TOXPROT()

            proteins
                .combine(DOWNLOAD_TOXPROT.out.blastdb)
                .set { toxprot_blastp_inputs }

            TOXPROT_BLASTP(toxprot_blastp_inputs)
        }

    emit:
        gff          = TRANSFER_FUNCTIONAL_ANNOTATIONS.out
        blastp       = BLASTP.out.tsv
        interproscan = INTERPROSCAN.out.tsv
        eggnog       = EGGNOG_MAPPER.out.annotations

        signalp = params.signalp
            ? SIGNALP.out.summary
            : Channel.empty()

        toxprot = params.toxprot_blast
            ? TOXPROT_BLASTP.out
            : Channel.empty()
}