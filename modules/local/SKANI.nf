
process SKANI{
	container = "quay.io/biocontainers/skani:0.0.1--h9f5acd7_0"

    input:
    path(fasta)
    path(skani_db)

    output:
    path("all.skani.tsv"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    """
    skani search ${fasta} -d ${skani_db} > all.skani.tsv
    """
} 
