include { GENERATE_RUNLIST    } from "$projectDir/modules/genes_annotation/rnaseq_assembly/generate_runlist.nf"
include { DOWNLOAD_SRA        } from "$projectDir/modules/genes_annotation/rnaseq_assembly/download_sra.nf"
include { STAR_SINGLE         } from "$projectDir/modules/genes_annotation/rnaseq_assembly/star_single.nf"
include { STAR_PAIRED         } from "$projectDir/modules/genes_annotation/rnaseq_assembly/star_paired.nf"
include { MERGE_BAM           } from "$projectDir/modules/genes_annotation/rnaseq_assembly/merge_bam.nf"
include { VARUS               } from "$projectDir/modules/genes_annotation/rnaseq_assembly/varus.nf"
include { SAMTOOLS_READ_COUNT } from "$projectDir/modules/genes_annotation/rnaseq_assembly/samtools_read_count.nf"


workflow RNASEQ_ASSEMBLY {
    take:
        genomes

    main:
        // Generate one SRA runlist per genome.
        genomes
            .map { id, meta, genome -> [id, meta] }
            .set { genomes_without_fasta }

        GENERATE_RUNLIST(genomes_without_fasta)

        // Keep filtered runlists containing at least one SRA accession.
        GENERATE_RUNLIST.out.star
            .filter { id, meta, runlist ->
                runlist.readLines().size() > 1
            }
            .set { star_runlists }

        // Keep complete runlists containing at least one SRA accession.
        GENERATE_RUNLIST.out.varus
            .filter { id, meta, runlist ->
                runlist.readLines().size() > 1
            }
            .set { varus_runlists }

        // Extract SRA accessions selected for STAR.
        star_runlists
            .flatMap { id, meta, runlist ->
                runlist.readLines()
                    .drop(1)
                    .findAll { it.trim() }
                    .collect { line ->
                        [id, line.split(/\s+/)[0]]
                    }
            }
            .set { sra_ids }

        if (params.rna_seq_aligner == 'star' || params.rna_seq_aligner == 'both') {
            log.info 'Running STAR for RNA-seq alignment'

            DOWNLOAD_SRA(sra_ids)

            // Extract local single-end reads.
            genomes
                .flatMap { id, meta, genome ->
                    def singleFiles = meta['RNASeq-file']?.single_end ?: []

                    singleFiles.collect { fastq ->
                        [id, fastq]
                    }
                }
                .set { rnaseq_single_files }

            // Extract local paired-end reads.
            genomes
                .flatMap { id, meta, genome ->
                    def pairedFiles = meta['RNASeq-file']?.paired_end ?: []

                    pairedFiles.collect { pair ->
                        [id, pair]
                    }
                }
                .set { rnaseq_paired_files }

            // Separate downloaded single-end and paired-end reads.
            DOWNLOAD_SRA.out
                .branch { id, sraFiles ->
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
                .map { id, fastqs, meta, genome ->
                    [id, meta, genome, fastqs.flatten()]
                }
                .set { star_single_inputs }

            // Combine downloaded and local paired-end reads.
            downloaded_sra_inputs.paired
                .concat(rnaseq_paired_files)
                .groupTuple(by: 0)
                .join(genomes, by: 0)
                .map { id, fastqs, meta, genome ->
                    [id, meta, genome, fastqs.flatten()]
                }
                .set { star_paired_inputs }

            STAR_SINGLE(star_single_inputs)
            STAR_PAIRED(star_paired_inputs)

            // Merge STAR BAM files when several read sets are available.
            STAR_SINGLE.out.bam
                .concat(STAR_PAIRED.out.bam)
                .groupTuple(by: [0, 1])
                .branch { id, meta, bams ->
                    multiple: bams.size() > 1
                    single: bams.size() == 1
                }
                .set { star_bam_groups }

            MERGE_BAM(star_bam_groups.multiple)

            star_bam_groups.single
                .map { id, meta, bams ->
                    [id, meta, bams.first()]
                }
                .set { single_star_bams }

            MERGE_BAM.out
                .concat(single_star_bams)
                .map { id, meta, bam ->
                    [id, meta, 'star', bam]
                }
                .set { all_star_bams }
        } else {
            all_star_bams = Channel.empty()
        }

        if (params.rna_seq_aligner == 'varus' || params.rna_seq_aligner == 'both') {
            log.info 'Running VARUS for RNA-seq alignment'

            // Associate each genome with its SRA runlist.
            genomes
                .join(varus_runlists, by: [0, 1])
                .set { genomes_for_varus }

            VARUS(genomes_for_varus)

            VARUS.out.bam
                .map { id, meta, bam ->
                    [id, meta, 'varus', bam]
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
            .branch { id, meta, aligner, readCountFile ->
                def readCount = readCountFile.readLines().first().toInteger()

                good_quality: readCount > 1000000
                poor_quality: readCount <= 1000000
            }
            .set { read_count_results }

        // Recover BAM files passing the read-count threshold.
        read_count_results.good_quality
            .join(all_bams, by: [0, 1, 2])
            .map { id, meta, aligner, readCountFile, bam ->
                [id, meta, aligner, bam]
            }
            .set { bams_for_gene_annotation }

    emit:
        bams = bams_for_gene_annotation
}