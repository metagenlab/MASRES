//Formatting fasta headers with AWK

process ASSEMBLY_HEADER_FORMAT {
	container = "docker://metagenlab/diag-pipeline-python-r:1.1"

	tag "$meta.id"

	publishDir "${params.outdir}/${meta.id}/04_hybrid_assembly", mode: 'copy'

	input:
	tuple val(meta), file(assembly)

	output:
	tuple val(meta), path("./${meta.id}_final_assembly.fasta"), emit: formatted_assembly
	path("./${meta.id}_final_assembly.fasta"), emit: ref_assembly
	shell:
	'''
	awk '/^>/{print ">contig_" ++i; next}{print}' < !{assembly} > !{meta.id}_final_assembly.fasta
	'''
}
