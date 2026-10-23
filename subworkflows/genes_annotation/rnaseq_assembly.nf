include { SUBSAMPLE_SRA_READS } from "$projectDir/modules/local/genes_annotation/rnaseq_assembly/subsample_sra_reads.nf"
include { GENERATE_RUNLIST    } from "$projectDir/modules/local/genes_annotation/rnaseq_assembly/generate_runlist.nf"
include { CUSTOM_SRATOOLSNCBISETTINGS } from "$projectDir/modules/nf-core/custom/sratoolsncbisettings/main.nf"
include { SRATOOLS_PREFETCH } from "$projectDir/modules/nf-core/sratools/prefetch/main.nf"
include { SRATOOLS_FASTERQDUMP } from "$projectDir/modules/nf-core/sratools/fasterqdump/main.nf"
include { STAR_GENOMEGENERATE } from "$projectDir/modules/nf-core/star/genomegenerate/main.nf"
include { STAR_ALIGN as STAR_SINGLE } from "$projectDir/modules/nf-core/star/align/main.nf"
include { STAR_ALIGN as STAR_PAIRED } from "$projectDir/modules/nf-core/star/align/main.nf"
include { SAMTOOLS_MERGE } from "$projectDir/modules/nf-core/samtools/merge/main.nf"
include { VARUS               } from "$projectDir/modules/local/genes_annotation/rnaseq_assembly/varus.nf"
include { SAMTOOLS_READ_COUNT } from "$projectDir/modules/local/genes_annotation/rnaseq_assembly/samtools_read_count.nf"


workflow RNASEQ_ASSEMBLY {
    take:
        genomes

    main:
        def imported = LegacyCacheUtils.load(params.legacy_cache_manifest, params.cacheDir)
        // Generate one SRA runlist per genome.
        genomes
            .map { meta, genome -> meta }
            .set { genomes_without_fasta }

        GENERATE_RUNLIST(genomes_without_fasta)

        // Keep only aligners whose raw prediction was not imported.
        GENERATE_RUNLIST.out.star
            .filter { meta, runlist -> runlist.readLines().size() > 1 && !LegacyCacheUtils.hasPrediction(imported, meta.id, 'braker-with-star') }
            .set { star_runlists }
        GENERATE_RUNLIST.out.varus
            .filter { meta, runlist -> runlist.readLines().size() > 1 && !LegacyCacheUtils.hasPrediction(imported, meta.id, 'braker-with-varus') }
            .set { varus_runlists }

        // Extract SRA accessions selected for STAR.
        star_runlists
            .flatMap { meta, runlist ->
                runlist.readLines()
                    .drop(1)
                    .findAll { it.trim() }
                    .collect { line ->
                        def fields = line.split(/\s+/)
                        [meta, fields[0], fields[4] == '0']
                    }
            }
            .set { sra_ids }

        if (params.rna_seq_aligner == 'star' || params.rna_seq_aligner == 'both') {
            log.info 'Running STAR for RNA-seq alignment'

            // Store run metadata separately: download tasks need unique run IDs,
            // while downstream evidence remains associated with its genome.
            sra_run_meta = sra_ids.map { meta, accession, singleEnd -> [accession, meta] }
            unique_sra_ids = sra_ids.map { meta, accession, singleEnd -> [[id: accession, single_end: singleEnd], accession] }.unique { it[1] }
            CUSTOM_SRATOOLSNCBISETTINGS(unique_sra_ids.map { meta, accession -> accession }.collect())
            SRATOOLS_PREFETCH(unique_sra_ids, CUSTOM_SRATOOLSNCBISETTINGS.out.ncbi_settings, [])
            SRATOOLS_FASTERQDUMP(SRATOOLS_PREFETCH.out.sra, CUSTOM_SRATOOLSNCBISETTINGS.out.ncbi_settings, [])
            run_reads = SRATOOLS_FASTERQDUMP.out.reads.map { meta, reads -> [meta.id, reads] }
                .combine(sra_run_meta, by: 0)
                .map { accession, reads, meta -> [meta, reads] }
            if (params.sra_max_spots != null) {
                SUBSAMPLE_SRA_READS(run_reads)
                downloaded_reads = SUBSAMPLE_SRA_READS.out.reads
            } else {
                downloaded_reads = run_reads
            }

            // Extract local single-end reads.
            genomes
                .flatMap { meta, genome ->
                    def singleFiles = LegacyCacheUtils.hasPrediction(imported, meta.id, 'braker-with-star') ? [] : (meta['RNASeq-file']?.single_end ?: [])

                    singleFiles.collect { fastq ->
                        [meta, fastq]
                    }
                }
                .set { rnaseq_single_files }

            // Extract local paired-end reads.
            genomes
                .flatMap { meta, genome ->
                    def pairedFiles = LegacyCacheUtils.hasPrediction(imported, meta.id, 'braker-with-star') ? [] : (meta['RNASeq-file']?.paired_end ?: [])

                    pairedFiles.collect { pair ->
                        [meta, pair]
                    }
                }
                .set { rnaseq_paired_files }

            // Separate downloaded single-end and paired-end reads.
            downloaded_reads
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

            needed_genomes = star_single_inputs.concat(star_paired_inputs)
                .map { meta, genome, reads -> [meta, genome] }.unique { it[0].id }
            STAR_GENOMEGENERATE(needed_genomes, Channel.value([[:], []]))
            star_index = STAR_GENOMEGENERATE.out.index.map { meta, index -> [meta.id, index] }
            single_reads = star_single_inputs.map { meta, genome, reads -> [meta + [single_end: true], reads] }
            paired_reads = star_paired_inputs.map { meta, genome, reads -> [meta + [single_end: false], reads] }
            single_with_index = single_reads.map { meta, reads -> [meta.id, meta, reads] }.join(star_index, failOnDuplicate: true, failOnMismatch: false)
            paired_with_index = paired_reads.map { meta, reads -> [meta.id, meta, reads] }.join(star_index, failOnDuplicate: true, failOnMismatch: false)
            // multiMap broadcasts paired read/index channels in the same order.
            single_channels = single_with_index.multiMap { id, meta, reads, index ->
                reads: [meta, reads]
                index: [[id: id], index]
            }
            paired_channels = paired_with_index.multiMap { id, meta, reads, index ->
                reads: [meta, reads]
                index: [[id: id], index]
            }
            STAR_SINGLE(single_channels.reads, single_channels.index, Channel.value([[:], []]), true)
            STAR_PAIRED(paired_channels.reads, paired_channels.index, Channel.value([[:], []]), true)

            // Merge STAR BAM files when several read sets are available.
            STAR_SINGLE.out.bam_sorted_aligned
                .concat(STAR_PAIRED.out.bam_sorted_aligned)
                .map { meta, bam -> [meta.findAll { key, value -> key != 'single_end' }, bam] }
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

        // Keep BAM files passing the configured mapped-read threshold.
        SAMTOOLS_READ_COUNT.out
            .branch { meta, readCountFile ->
                def readCount = readCountFile.readLines().first().toInteger()

                good_quality: readCount > params.rnaseq_min_mapped_reads
                poor_quality: readCount <= params.rnaseq_min_mapped_reads
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
    versions = params.rna_seq_aligner in ['star', 'both'] ? SAMTOOLS_MERGE.out.versions.mix(STAR_GENOMEGENERATE.out.versions).mix(STAR_SINGLE.out.versions).mix(STAR_PAIRED.out.versions).mix(SRATOOLS_PREFETCH.out.versions).mix(SRATOOLS_FASTERQDUMP.out.versions).mix(CUSTOM_SRATOOLSNCBISETTINGS.out.versions) : Channel.empty()
        bams = bams_for_gene_annotation
}