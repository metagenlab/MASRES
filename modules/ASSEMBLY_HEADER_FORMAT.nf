//Formatting fasta headers with AWK

process ASSEMBLY_HEADER_FORMAT {
	container = "docker://metagenlab/diag-pipeline-python-r:1.1"

	tag "$ID"

	publishDir "${params.outdir}/${ID}/04_hybrid_assembly", mode: 'copy'

	input:
	tuple val(ID), file(assembly)

	output:
	val ID, emit: ID
	tuple val(ID), path("./${ID}_final_assembly.fasta"), emit: formatted_assembly

	shell:
	'''
	awk 'BEGIN{FS="_"}{if(/^>/){print $1"_"$2}else{print $0}}' !{assembly} > !{ID}_final_assembly.fasta
	'''
}
