// Assembly with flye

process FLYE {
        container = "docker://staphb/flye"

	tag "$meta.id"	

	publishDir "${params.outdir}/${meta.id}/", mode: 'copy'

        input:
        tuple val(meta), file(reads)

        output:
	val meta, emit: meta
	tuple val(meta), path(reads), path("./01_assembly_flye/assembly.fasta"), emit: assembly_channel
        tuple val(meta), path("./01_assembly_flye/assembly.fasta"), emit: assembly
	tuple val(meta), path("./01_assembly_flye/assembly_graph.gfa"), emit: graph_gfa
	tuple val(meta), path("./01_assembly_flye/assembly_info.txt"), emit: assembly_info

        script:
        
        """
	flye --nano-raw ${reads} --out-dir 01_assembly_flye --threads $params.threads
        """
}

