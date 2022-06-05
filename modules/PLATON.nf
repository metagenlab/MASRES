process PLATON {
        container = "docker://quay.io/biocontainers/platon:1.6--pyhdfd78af_1"

        tag "$ID"

        publishDir "${params.outdir}/${ID}/", mode: 'copy'

        input:

        tuple val(ID), file(assembly), file(depth)
        path db

        output:
        val ID, emit: ID
        tuple val(ID), file(depth), emit: depth
        tuple val(ID), path("./06_plasmid_annotation/${ID}.tsv"), emit: plasmid_annot

        shell:
        """
        platon --db ${db} --output 06_plasmid_annotation --prefix ${ID} ${assembly}
        """
}


