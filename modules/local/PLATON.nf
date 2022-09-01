process PLATON {
        container = "docker://quay.io/biocontainers/platon:1.6--pyhdfd78af_1"

        tag "$meta.id"

        publishDir "${params.outdir}/${meta.id}/", mode: 'copy'

        input:

        tuple val(meta), file(assembly)
        path db

        output:
        val meta, emit: meta
        tuple val(meta), path("./06_plasmid_annotation/plasmids-${meta.id}.tsv"), emit: plasmid_annot

        shell:
        """
        platon --db ${db} --output 06_plasmid_annotation --prefix plasmids-${meta.id} ${assembly}
        """
}


