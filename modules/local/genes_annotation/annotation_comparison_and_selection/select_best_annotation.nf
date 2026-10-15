process SELECT_BEST_ANNOTATION {
    label 'exogap_python'
    tag { 'thresholds_v1_' + java.security.MessageDigest.getInstance('SHA-256').digest(groovy.json.JsonOutput.toJson(new TreeMap(params.selection_thresholds)).getBytes('UTF-8')).encodeHex().toString().take(16) }
    input:
    tuple val(status), path(summary)
    output:
    tuple val(status), path('best_annotations.tsv'), emit: selection
    path 'selection_scores.tsv', emit: scores
    path 'selection_thresholds.json', emit: thresholds
    script:
    def thresholdsJson = groovy.json.JsonOutput.toJson(params.selection_thresholds)
    """
    cat > selection_config.json <<'SELECTION_CONFIG'
    ${thresholdsJson}
    SELECTION_CONFIG
    select_best_annotations.py -i ${summary} -o best_annotations.tsv \
        --thresholds selection_config.json --scores-output selection_scores.tsv \
        --thresholds-output selection_thresholds.json
    """
    stub:
    """
    touch best_annotations.tsv selection_scores.tsv
    echo '{}' > selection_thresholds.json
    """
}
