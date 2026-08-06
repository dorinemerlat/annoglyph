process PUBLISH_BEST_ANNOTATIONS {

    tag "$id"

    input:
    tuple val(id), val(meta), path(file)

    output:
    tuple val(id), val(meta), path("${id}.*")

    script:
    def ext = file.getExtension()

    """
    cp ${file} ${id}.${ext}
    """
}