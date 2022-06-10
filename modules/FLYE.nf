// Assembly with flye

process FLYE {
        container = "docker://staphb/flye"

	tag "$ID"	

	publishDir "${params.outdir}/${ID}/", mode: 'copy'

        input:
        tuple val(ID), file(reads)

        output:
	val ID, emit: ID
        tuple val(ID), path("./01_assembly_flye/assembly.fasta"), emit: assembly
	tuple val(ID), path("./01_assembly_flye/assembly_graph.gfa"), emit: graph_gfa
	tuple val(ID), path("./01_assembly_flye/assembly_info.txt"), emit: assembly_info

        script:
        
        """
	flye --nano-raw ${reads} --out-dir 01_assembly_flye --threads $params.threads
        """
}

