include { INTERPROSCAN } from './modules/nf-core/interproscan/main.nf'
workflow {
    INTERPROSCAN(Channel.of([[id:'ips_smoke', annotation_method:'braker2'], file("${projectDir}/test/validation/genes_real.fasta", checkIfExists: true)]), params.interproscan_db ? file(params.interproscan_db, checkIfExists: true) : [])
    INTERPROSCAN.out.tsv.view { meta, tsv ->
        assert tsv.exists()
        "INTERPROSCAN_REAL_OK: ${tsv.name} (${tsv.size()} bytes)"
    }
}
