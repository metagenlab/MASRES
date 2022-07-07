//Polish with Homopolish


process HOMOPOLISH {
        container = "docker://quay.io/biocontainers/homopolish:0.3.3--pyh5e36f6f_0"

        tag "$meta.id"

        publishDir "${params.outdir}/${meta.id}/03_polished", mode: 'copy'

        input:
        tuple val(meta), file(medaka_polish)
	path db

        output:
        val(meta), emit: meta
        tuple val(meta), file(polished_assembly), emit: homopolished

        script:
	polished_assembly="${meta.id}_homopolished.fasta"
        """
	mv ${medaka_polish} ${meta.id}.fasta
        homopolish polish -a ${meta.id}.fasta -s ${db} -o . -m R9.4.pkl -t ${params.threads}
	"""
}
