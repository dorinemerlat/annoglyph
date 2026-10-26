include { SUBSAMPLE_SRA_READS } from "$projectDir/modules/local/genes_annotation/rnaseq_assembly/subsample_sra_reads.nf"
include { GENERATE_RUNLIST    } from "$projectDir/modules/local/genes_annotation/rnaseq_assembly/generate_runlist.nf"
include { CUSTOM_SRATOOLSNCBISETTINGS } from "$projectDir/modules/nf-core/custom/sratoolsncbisettings/main.nf"
include { SRATOOLS_PREFETCH } from "$projectDir/modules/nf-core/sratools/prefetch/main.nf"
include { SRATOOLS_FASTERQDUMP } from "$projectDir/modules/nf-core/sratools/fasterqdump/main.nf"
include { STAR_GENOMEGENERATE } from "$projectDir/modules/nf-core/star/genomegenerate/main.nf"
include { STAR_ALIGN as STAR_SINGLE } from "$projectDir/modules/nf-core/star/align/main.nf"
include { STAR_ALIGN as STAR_PAIRED } from "$projectDir/modules/nf-core/star/align/main.nf"
include { CLEAN_RNASEQ_READS } from "$projectDir/modules/local/genes_annotation/rnaseq_assembly/clean_reads.nf"
include { SORT_RNASEQ_BAM } from "$projectDir/modules/local/genes_annotation/rnaseq_assembly/sort_bam.nf"
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
                .map { accession, reads, meta ->
                    def ordered = (reads instanceof List ? reads : [reads]).sort { it.name }
                    [RnaseqUtils.readset(meta, ordered, 'sra', accession), ordered]
                }
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
                        [RnaseqUtils.readset(meta, fastq, 'local_se'), fastq]
                    }
                }
                .set { rnaseq_single_files }

            // Extract local paired-end reads.
            genomes
                .flatMap { meta, genome ->
                    def pairedFiles = LegacyCacheUtils.hasPrediction(imported, meta.id, 'braker-with-star') ? [] : (meta['RNASeq-file']?.paired_end ?: [])

                    pairedFiles.collect { pair ->
                        [RnaseqUtils.readset(meta, pair, 'local_pe'), pair]
                    }
                }
                .set { rnaseq_paired_files }

            // Each SRA run or local library remains independent through QC,
            // alignment and sorting. The genome identity is restored at the merge.
            readsets = downloaded_reads.concat(rnaseq_single_files, rnaseq_paired_files)
                .unique { meta, reads -> meta.id }
            CLEAN_RNASEQ_READS(readsets)
            clean_reads = CLEAN_RNASEQ_READS.out.reads
                .map { meta, reads, report ->
                    def ordered = (reads instanceof List ? reads : [reads]).sort { it.name }
                    [meta + RnaseqUtils.readStats(report), ordered]
                }
                .filter { meta, reads ->
                    if (meta.rnaseq_read_count == 0) log.warn "No usable RNA-seq reads after QC: ${meta.id}; retaining its QC report"
                    meta.rnaseq_read_count > 0
                }
            star_qc_reports = CLEAN_RNASEQ_READS.out.reads.map { meta, reads, report -> [meta, report] }

            // combine broadcasts the single genome/index to every run; join
            // would consume a genome key once and lose the remaining libraries.
            read_genomes = clean_reads.map { meta, reads -> [meta.genome_id, meta, reads] }
                .combine(genomes.map { meta, genome -> [meta.id, genome] }, by: 0)
            needed_genomes = read_genomes.map { id, meta, reads, genome -> [meta.genome_meta, genome] }
                .unique { meta, genome -> meta.id }
            STAR_GENOMEGENERATE(needed_genomes, Channel.value([[:], []]))
            star_index = STAR_GENOMEGENERATE.out.index.map { meta, index -> [meta.id, index, RnaseqUtils.indexBytes(index)] }
            aligned_inputs = clean_reads.map { meta, reads -> [meta.genome_id, meta, reads] }
                .combine(star_index, by: 0)
                .map { id, meta, reads, index, bytes ->
                    [meta + [star_memory_gb: RnaseqUtils.starMemoryGb(bytes, meta.rnaseq_read_count as long, params)], reads, index]
                }
                .branch {
                    single: it[0].single_end
                    paired: !it[0].single_end
                }
            single_channels = aligned_inputs.single.multiMap { meta, reads, index ->
                reads: [meta, reads]
                index: [[id: meta.genome_id], index]
            }
            paired_channels = aligned_inputs.paired.multiMap { meta, reads, index ->
                reads: [meta, reads]
                index: [[id: meta.genome_id], index]
            }
            STAR_SINGLE(single_channels.reads, single_channels.index, Channel.value([[:], []]), true)
            STAR_PAIRED(paired_channels.reads, paired_channels.index, Channel.value([[:], []]), true)
            SORT_RNASEQ_BAM(STAR_SINGLE.out.bam.concat(STAR_PAIRED.out.bam))
            star_run_bams = SORT_RNASEQ_BAM.out.bam

            // Include the library inventory in the merge cache key and sort it
            // deterministically, independently of task completion order.
            star_bam_groups = SORT_RNASEQ_BAM.out.bam
                .map { meta, bam -> [meta.genome_id, meta, bam] }
                .groupTuple(by: 0)
                .map { id, metas, bams ->
                    def ordered = (0..<bams.size()).toList().sort { a, b -> metas[a].readset_id <=> metas[b].readset_id }
                    [RnaseqUtils.mergeMeta(metas), ordered.collect { bams[it] }]
                }
                .branch {
                    multiple: it[1].size() > 1
                    single: it[1].size() == 1
                }
            SAMTOOLS_MERGE(star_bam_groups.multiple.map { meta, bams -> [meta, bams, []] }, Channel.value([[:], [], [], []]), '')
            single_star_bams = star_bam_groups.single.map { meta, bams -> [meta, bams.first()] }
            all_star_bams = SAMTOOLS_MERGE.out.bam.concat(single_star_bams)
                .map { meta, bam -> [meta.genome_meta + [rnaseq_aligner: 'star'], bam] }
            star_versions = CLEAN_RNASEQ_READS.out.versions.mix(SORT_RNASEQ_BAM.out.versions,
                SAMTOOLS_MERGE.out.versions, STAR_GENOMEGENERATE.out.versions,
                STAR_SINGLE.out.versions, STAR_PAIRED.out.versions,
                SRATOOLS_PREFETCH.out.versions, SRATOOLS_FASTERQDUMP.out.versions,
                CUSTOM_SRATOOLSNCBISETTINGS.out.versions)
        } else {
            all_star_bams = Channel.empty()
            star_versions = Channel.empty()
            star_run_bams = Channel.empty()
            star_qc_reports = Channel.empty()
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
                def readCount = readCountFile.readLines().first().toLong()

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
        versions = star_versions
        bams = bams_for_gene_annotation
        run_bams = star_run_bams
        qc_reports = star_qc_reports
}
