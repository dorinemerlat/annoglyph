include { BLASTP                          } from "$projectDir/modules/local/genes_annotation/functional_annotation/blastp.nf"
include { INTERPROSCAN } from "$projectDir/modules/nf-core/interproscan/main.nf"
include { EGGNOGMAPPER } from "$projectDir/modules/nf-core/eggnogmapper/main.nf"
include { TRANSFER_FUNCTIONAL_ANNOTATIONS } from "$projectDir/modules/local/genes_annotation/functional_annotation/transfer_functional_annotations.nf"
include { DOWNLOAD_TOXPROT                } from "$projectDir/modules/local/genes_annotation/functional_annotation/download_toxprot.nf"
include { TOXPROT_BLASTP                  } from "$projectDir/modules/local/genes_annotation/functional_annotation/toxprot_blastp.nf"
include { SIGNALP                         } from "$projectDir/modules/local/genes_annotation/functional_annotation/signalp.nf"


workflow FUNCTIONAL_ANNOTATION {

    take:
        gff
        proteins

    main:
        if (!params.eggnog_db) error 'Configure params.eggnog_db with a versioned eggNOG database directory in your site config.'
        // Search canonical proteins against Swiss-Prot.
        proteins
            .map { meta, fasta ->
                [meta, fasta, params.blast_db ?: 'uniprot_sprot', 1e-6, 3, 20]
            }
            .set { blastp_inputs }

        BLASTP(blastp_inputs)

        // Identify protein domains and functional signatures.
        interproscan_data = params.interproscan_db ? Channel.value(file(params.interproscan_db, checkIfExists: true)) : []
        INTERPROSCAN(proteins, interproscan_data)

        // Assign orthology-based functional annotations.
        EGGNOGMAPPER(proteins, Channel.value(['diamond', []]), Channel.value(file(params.eggnog_db, checkIfExists: true)))

        // Transfer functional annotations to the GFF3.
        gff
            .join(BLASTP.out.tsv, by: 0)
            .join(INTERPROSCAN.out.tsv, by: 0)
            .join(EGGNOGMAPPER.out.annotations, by: 0)
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
    versions = INTERPROSCAN.out.versions.mix(EGGNOGMAPPER.out.versions)
        gff          = TRANSFER_FUNCTIONAL_ANNOTATIONS.out
        blastp       = BLASTP.out.tsv
        interproscan = INTERPROSCAN.out.tsv
        eggnog       = EGGNOGMAPPER.out.annotations

        signalp = params.signalp
            ? SIGNALP.out.summary
            : Channel.empty()

        toxprot = params.toxprot_blast
            ? TOXPROT_BLASTP.out
            : Channel.empty()
}