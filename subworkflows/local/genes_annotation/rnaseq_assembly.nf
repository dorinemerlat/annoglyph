// Include local modules
include { GENERATE_RUNLIST      } from "$projectDir/modules/local/module_genes/rnaseq_assembly/generate_runlist.nf"
include { DOWNLOAD_SRA          } from "$projectDir/modules/local/module_genes/rnaseq_assembly/download_sra.nf"
include { STAR_SINGLE           } from "$projectDir/modules/local/module_genes/rnaseq_assembly/star_single.nf"
include { STAR_PAIRED           } from "$projectDir/modules/local/module_genes/rnaseq_assembly/star_paired.nf"
include { MERGE_BAM             } from "$projectDir/modules/local/module_genes/rnaseq_assembly/merge_bam.nf"
include { VARUS                 } from "$projectDir/modules/local/module_genes/rnaseq_assembly/varus.nf"
include { SAMTOOLS_READ_COUNT   } from "$projectDir/modules/local/module_genes/rnaseq_assembly/samtools_read_count.nf"


workflow RNASEQ_ASSEMBLY {

    take:
        /*
         * Input genomes from the preprocessing step.
         * These are unmasked genomes and are used for STAR or VARUS RNA-seq alignment.
         *
         * Expected format:
         *   [id, meta, fasta]
         */
        genomes

    main:

        /*
         * Generate one SRA runlist per genome/species.
         *
         * The genome FASTA is not needed at this step, so only the id and metadata
         * are passed to GENERATE_RUNLIST.
         */
        genomes
            .map { id, meta, fasta -> [ id, meta ] }
            .set { genomes_without_fasta }

        GENERATE_RUNLIST(genomes_without_fasta)


        /*
         * Keep only non-empty runlists.
         *
         * A valid runlist contains a header line plus at least one SRA run.
         */
        GENERATE_RUNLIST.out
            .filter { id, meta, runlist ->
                runlist.readLines().size() > 1
            }
            .set { runlists }


        /*
         * Extract individual SRA run accessions from each runlist.
         *
         * Output format:
         *   [id, sra_id]
         */
        runlists
            .flatMap { id, meta, runlist ->
                runlist.readLines()
                    .drop(1)
                    .findAll { it.trim() }
                    .collect { line ->
                        [ id, line.split(/\s+/)[0] ]
                    }
            }
            .set { sra_IDs }


        /*
         * STAR alignment branch.
         *
         * This branch combines:
         *   - SRA-derived FASTQ files downloaded by DOWNLOAD_SRA
         *   - user-provided local RNA-seq FASTQ files from meta['RNASeq-file']
         *
         * STAR is run separately for single-end and paired-end reads.
         */
        if (params.rna_seq_aligner == 'star' || params.rna_seq_aligner == 'both') {

            log.info "Running STAR for RNA-seq alignment"

            /*
             * Download SRA FASTQ files for all selected SRA accessions.
             */
            DOWNLOAD_SRA(sra_IDs)


            /*
             * Extract local single-end RNA-seq FASTQ files from the metadata.
             *
             * Expected metadata structure:
             *   meta['RNASeq-file'].single_end = [file1, file2, ...]
             *
             * Output format:
             *   [id, fastq]
             */
            genomes
                .flatMap { id, meta, genome ->
                    def rnaseq = meta['RNASeq-file']
                    def single_files = rnaseq?.single_end ?: []

                    single_files.collect { fastq ->
                        [ id, fastq ]
                    }
                }
                .set { rnaseq_single_files }


            /*
             * Extract local paired-end RNA-seq FASTQ files from the metadata.
             *
             * Expected metadata structure:
             *   meta['RNASeq-file'].paired_end = [[r1, r2], [r1, r2], ...]
             *
             * Output format:
             *   [id, [r1, r2]]
             */
            genomes
                .flatMap { id, meta, genome ->
                    def rnaseq = meta['RNASeq-file']
                    def paired_files = rnaseq?.paired_end ?: []

                    paired_files.collect { pair ->
                        [ id, pair ]
                    }
                }
                .set { rnaseq_paired_files }


            /*
             * Split downloaded SRA FASTQ files according to read layout.
             *
             * DOWNLOAD_SRA emits either:
             *   [id, fastq]
             * or:
             *   [id, [r1, r2]]
             */
            DOWNLOAD_SRA.out
                .branch { id, sra_files ->
                    def files = (sra_files instanceof List) ? sra_files : [sra_files]

                    single: files.size() == 1
                    paired: files.size() == 2
                }
                .set { downloaded_sra_input }


            /*
             * Build STAR input for single-end reads by combining downloaded SRA
             * reads and local user-provided reads.
             *
             * Output format:
             *   [id, meta, genome, fastq_list]
             */
            downloaded_sra_input.single
                .concat(rnaseq_single_files)
                .groupTuple(by: 0)
                .join(genomes, by: 0)
                .map { id, fastq, meta, genome ->
                    [ id, meta, genome, fastq.flatten() ]
                }
                .set { star_input_single_end }


            /*
             * Build STAR input for paired-end reads by combining downloaded SRA
             * reads and local user-provided reads.
             *
             * Output format:
             *   [id, meta, genome, fastq_pair_list]
             */
            downloaded_sra_input.paired
                .concat(rnaseq_paired_files)
                .groupTuple(by: 0)
                .join(genomes, by: 0)
                .map { id, fastq, meta, genome ->
                    [ id, meta, genome, fastq.flatten() ]
                }
                .set { star_input_paired_end }


            /*
             * Run STAR on single-end and paired-end RNA-seq inputs.
             */
            STAR_SINGLE(star_input_single_end)
            STAR_PAIRED(star_input_paired_end)


            /*
             * Combine all STAR BAM files per genome.
             *
             * If several BAM files are produced for the same genome, merge them.
             * If only one BAM file exists, keep it directly.
             */
            STAR_SINGLE.out.bam
                .concat(STAR_PAIRED.out.bam)
                .groupTuple(by: [0, 1])
                .map { id, meta, bams ->
                    [ id, meta, bams ]
                }
                .branch { id, meta, bams ->
                    has_multiple_bams: bams.size() > 1
                    single_bam: bams.size() == 1
                }
                .set { merged_star_input }

            MERGE_BAM(merged_star_input.has_multiple_bams)

            MERGE_BAM.out
                .concat(merged_star_input.single_bam)
                .map { id, meta, bam ->
                    [ id, meta, 'star', bam ]
                }
                .set { all_star_bams }

        }
        else {
            /*
             * Empty channel used when STAR is disabled.
             */
            Channel.empty().set { all_star_bams }
        }


        /*
         * VARUS alignment branch.
         *
         * VARUS uses the genome and SRA runlist directly.
         */
        if (params.rna_seq_aligner == 'varus' || params.rna_seq_aligner == 'both') {

            log.info "Running VARUS for RNA-seq alignment"

            /*
             * Join genomes with their non-empty runlists.
             *
             * Output format:
             *   [id, meta, genome, runlist]
             */
            genomes
                .join(runlists, by: [0, 1])
                .set { genomes_for_varus }

            VARUS(genomes_for_varus)

            VARUS.out.bam
                .map { id, meta, bam ->
                    [ id, meta, 'varus', bam ]
                }
                .set { all_varus_bams }

        }
        else {
            /*
             * Empty channel used when VARUS is disabled.
             */
            Channel.empty().set { all_varus_bams }
        }


        /*
         * Combine BAM files generated by all enabled aligners.
         *
         * Output format:
         *   [id, meta, aligner, bam]
         */
        all_star_bams
            .concat(all_varus_bams)
            .set { all_bams }


        /*
         * Count mapped reads in each BAM file.
         */
        SAMTOOLS_READ_COUNT(all_bams)


        /*
         * Keep only BAM files with enough mapped reads for gene annotation.
         *
         * Current threshold:
         *   > 1,000,000 mapped reads
         */
        SAMTOOLS_READ_COUNT.out
            .branch { id, meta, aligner, read_count_file ->
                def read_number = read_count_file.readLines()[0].toInteger()

                good_quality: read_number > 1000000
                poor_quality: read_number <= 1000000
            }
            .set { read_count_results }


        /*
         * Recover the corresponding BAM files for the good-quality read counts.
         *
         * Output format:
         *   [id, meta, aligner, bam]
         */
        read_count_results.good_quality
            .join(all_bams, by: [0, 1, 2])
            .map { id, meta, aligner, read_count_file, bam ->
                [ id, meta, aligner, bam ]
            }
            .set { bams_for_gene_annotation }


    emit:
        /*
         * BAM files suitable for gene annotation.
         *
         * Output format:
         *   [id, meta, aligner, bam]
         */
        bams = bams_for_gene_annotation
}