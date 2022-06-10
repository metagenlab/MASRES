//Polish with Homopolish


process HOMOPOLISH {
        container = "docker://quay.io/biocontainers/homopolish:0.3.3--pyh5e36f6f_0"

        tag "$ID"

        publishDir "${params.outdir}/${ID}/", mode: 'copy'

        input:
        tuple val(ID), file(medaka_polish)
	path db

        output:
        val(ID), emit: ID
        tuple val(ID), path('./03_homopolish/consensus_homopolished.fasta'), emit: homopolished

        script:
        """
        homopolish polish -a ${medaka_polish} -s ${db} -o 03_homopolish -m R9.4.pkl -t ${params.threads}
        """
}
