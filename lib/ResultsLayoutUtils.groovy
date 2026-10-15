class ResultsLayoutUtils {
    static List assets(meta, files, directory) {
        def paths = files instanceof List ? files.flatten() : [files]
        paths.collect { path -> [meta.id, [path, "${directory}/${path.name}"]] }
    }
}
