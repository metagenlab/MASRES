process SOURMASH_GATHER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::4.4.3--hdfd78af_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.4.3--hdfd78af_0':
        'quay.io/biocontainers/sourmash:4.4.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(signature)
    path(sketch_db)

    output:
    tuple val(meta), path("*.mash"), emit: csv_tuple
    path("*.mash")                , emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sourmash gather \\
        $args \\
	-o '${prefix}.mash' \\
        ${signature} \\
        ${sketch_db}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}
