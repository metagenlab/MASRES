process QUAST {
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::quast=5.2.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'quay.io/biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"
    
    tag "${meta.id}"


    input:
    tuple val(meta), path(assembly)
    path gff
    val use_fasta
    val use_gff

    output:
    path "*"    , emit: results
    tuple val(meta), path("./report.tsv"), emit: tsv_tuple
    path "*.tsv"        , emit: tsv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    quast.py \\
        --output-dir . \\
        --threads $task.cpus \\
        $args \\
        ${assembly}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
