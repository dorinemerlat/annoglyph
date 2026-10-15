include { GENERATE_RUNLIST    } from "$projectDir/modules/local/genes_annotation/rnaseq_assembly/generate_runlist.nf"
include { DOWNLOAD_SRA        } from "$projectDir/modules/local/genes_annotation/rnaseq_assembly/download_sra.nf"
include { STAR_SINGLE         } from "$projectDir/modules/local/genes_annotation/rnaseq_assembly/star_single.nf"
include { STAR_PAIRED         } from "$projectDir/modules/local/genes_annotation/rnaseq_assembly/star_paired.nf"
include { SAMTOOLS_MERGE } from "$projectDir/modules/nf-core/samtools/merge/main.nf"
include { VARUS               } from "$projectDir/modules/local/genes_annotation/rnaseq_assembly/varus.nf"
include { SAMTOOLS_READ_COUNT } from "$projectDir/modules/local/genes_annotation/rnaseq_assembly/samtools_read_count.nf"


workflow RNASEQ_ASSEMBLY {
    take:
        genomes

    main:
        // Generate one SRA runlist per genome.
        genomes
            .map { meta, genome -> meta }
            .set { genomes_without_fasta }

        GENERATE_RUNLIST(genomes_without_fasta)

        // Keep filtered runlists containing at least one SRA accession.
        GENERATE_RUNLIST.out.star
            .filter { meta, runlist ->
                runlist.readLines().size() > 1
            }
            .set { star_runlists }

        // Keep complete runlists containing at least one SRA accession.
        GENERATE_RUNLIST.out.varus
            .filter { meta, runlist ->
                runlist.readLines().size() > 1
            }
            .set { varus_runlists }

        // Extract SRA accessions selected for STAR.
        star_runlists
            .flatMap { meta, runlist ->
                runlist.readLines()
                    .drop(1)
                    .findAll { it.trim() }
                    .collect { line ->
                        [meta, line.split(/\s+/)[0]]
                    }
            }
            .set { sra_ids }

        if (params.rna_seq_aligner == 'star' || params.rna_seq_aligner == 'both') {
            log.info 'Running STAR for RNA-seq alignment'

            DOWNLOAD_SRA(sra_ids)

            // Extract local single-end reads.
            genomes
                .flatMap { meta, genome ->
                    def singleFiles = meta['RNASeq-file']?.single_end ?: []

                    singleFiles.collect { fastq ->
                        [meta, fastq]
                    }
                }
                .set { rnaseq_single_files }

            // Extract local paired-end reads.
            genomes
                .flatMap { meta, genome ->
                    def pairedFiles = meta['RNASeq-file']?.paired_end ?: []

                    pairedFiles.collect { pair ->
                        [meta, pair]
                    }
                }
                .set { rnaseq_paired_files }

            // Separate downloaded single-end and paired-end reads.
            DOWNLOAD_SRA.out
                .branch { meta, sraFiles ->
                    def files = sraFiles instanceof List ? sraFiles : [sraFiles]

                    single: files.size() == 1
                    paired: files.size() == 2
                }
                .set { downloaded_sra_inputs }

            // Combine downloaded and local single-end reads.
            downloaded_sra_inputs.single
                .concat(rnaseq_single_files)
                .groupTuple(by: 0)
                .join(genomes, by: 0)
                .map { meta, fastqs, genome ->
                    [meta, genome, fastqs.flatten()]
                }
                .set { star_single_inputs }

            // Combine downloaded and local paired-end reads.
            downloaded_sra_inputs.paired
                .concat(rnaseq_paired_files)
                .groupTuple(by: 0)
                .join(genomes, by: 0)
                .map { meta, fastqs, genome ->
                    [meta, genome, fastqs.flatten()]
                }
                .set { star_paired_inputs }

            STAR_SINGLE(star_single_inputs)
            STAR_PAIRED(star_paired_inputs)

            // Merge STAR BAM files when several read sets are available.
            STAR_SINGLE.out.bam
                .concat(STAR_PAIRED.out.bam)
                .groupTuple(by: 0)
                .branch { meta, bams ->
                    multiple: bams.size() > 1
                    single: bams.size() == 1
                }
                .set { star_bam_groups }

            SAMTOOLS_MERGE(star_bam_groups.multiple.map { meta, bams -> [meta, bams, []] }, Channel.value([[:], [], [], []]), '')

            star_bam_groups.single
                .map { meta, bams ->
                    [meta, bams.first()]
                }
                .set { single_star_bams }

            SAMTOOLS_MERGE.out.bam
                .concat(single_star_bams)
                .map { meta, bam ->
                    [meta + [rnaseq_aligner: 'star'], bam]
                }
                .set { all_star_bams }
        } else {
            all_star_bams = Channel.empty()
        }

        if (params.rna_seq_aligner == 'varus' || params.rna_seq_aligner == 'both') {
            log.info 'Running VARUS for RNA-seq alignment'

            // Associate each genome with its SRA runlist.
            genomes
                .join(varus_runlists, by: 0)
                .set { genomes_for_varus }

            VARUS(genomes_for_varus)

            VARUS.out.bam
                .map { meta, bam ->
                    [meta + [rnaseq_aligner: 'varus'], bam]
                }
                .set { all_varus_bams }
        } else {
            all_varus_bams = Channel.empty()
        }

        // Combine BAM files from all enabled aligners.
        all_star_bams
            .concat(all_varus_bams)
            .set { all_bams }

        SAMTOOLS_READ_COUNT(all_bams)

        // Keep BAM files with more than one million mapped reads.
        SAMTOOLS_READ_COUNT.out
            .branch { meta, readCountFile ->
                def readCount = readCountFile.readLines().first().toInteger()

                good_quality: readCount > 1000000
                poor_quality: readCount <= 1000000
            }
            .set { read_count_results }

        // Recover BAM files passing the read-count threshold.
        read_count_results.good_quality
            .join(all_bams, by: 0)
            .map { meta, readCountFile, bam ->
                [meta, bam]
            }
            .set { bams_for_gene_annotation }

    emit:
    versions = params.rna_seq_aligner in ['star', 'both'] ? SAMTOOLS_MERGE.out.versions : Channel.empty()
        bams = bams_for_gene_annotation
}