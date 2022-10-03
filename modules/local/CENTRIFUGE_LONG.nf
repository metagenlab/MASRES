process CENTRIFUGE_LONG {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::centrifuge=1.0.4_beta" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/centrifuge:1.0.4_beta--h9a82719_6' :
        'quay.io/biocontainers/centrifuge:1.0.4_beta--h9a82719_6' }"

    input:
    tuple val(meta), path(reads)
    path db
    val save_unaligned
    val save_aligned
    val sam_format

    output:
    tuple val(meta), path('*report.txt')                 , emit: report
    tuple val(meta), path('*results.txt')                , emit: results
    tuple val(meta), path('*.sam')                       , optional: true, emit: sam
    tuple val(meta), path('*.mapped.fastq{,.1,.2}.gz')   , optional: true, emit: fastq_mapped
    tuple val(meta), path('*.unmapped.fastq{,.1,.2}.gz') , optional: true, emit: fastq_unmapped
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = "-U ${reads}"
    def unaligned = ''
    def aligned = ''
    
    unaligned = save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ''
    aligned = save_aligned ? "--al-gz ${prefix}.mapped.fastq.gz" : ''

    def sam_output = sam_format ? "--out-fmt 'sam'" : ''
    """
    ## we add "-no-name ._" to ensure silly Mac OSX metafiles files aren't included
    db_name=`find -L ${db} -name "*.1.cf" -not -name "._*"  | sed 's/.1.cf//'`
    centrifuge \\
        -x \$db_name \\
        -p $task.cpus \\
        $paired \\
        --report-file ${prefix}.report.txt \\
        -S ${prefix}_long.results.txt \\
        $unaligned \\
        $aligned \\
        $sam_output \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuge: \$( centrifuge --version  | sed -n 1p | sed 's/^.*centrifuge-class version //')
    END_VERSIONS
    """
}
