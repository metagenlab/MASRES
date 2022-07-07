//Polish with Medaka


process MEDAKA {
        container = "docker://ontresearch/medaka"

	tag "$ID"
	
	publishDir "${params.outdir}/${ID}/", mode: 'copy'
	
        input:
        tuple val(ID), file(reads), file(assembly)

        output:
	val(ID), emit: ID
        tuple val(ID), path("./02_polishing_medaka/consensus.fasta"), emit: medaka_polish 

        script:
        """
        medaka_consensus -i ${reads} -d ${assembly} -o 02_polishing_medaka -t ${params.threads}
        """
}
