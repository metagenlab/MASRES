//Assembly with spades

process SPADES {

	container = "docker://staphb/spades"
	
	tag "$ID"

	publishDir "${params.outdir}/${ID}/", mode: 'copy'

	input:
	tuple val(ID), file(read1), file(read2)
	
	output:
	val ID, emit: ID
	tuple val(ID), path('./01_spades/scaffolds.fasta'), emit: scaffolds
	tuple val(ID), path('./01_spades/contigs.fasta'), emit: contigs
	tuple val(ID), path('./01_spades/assembly_graph_with_scaffolds.gfa'), emit: assembly_scaff
	tuple val(ID), path('./01_spades/assembly_graph.fastg'), emit: assembly_fastg

	script:
	"""
	set +u
	spades.py -k 21,33,55,77,99,127 -1 ${read1} -2 ${read2} -o 01_spades --isolate --threads ${params.threads}
	set -u
	"""
}
	
