//Formatting fasta headers with AWK

process ASSEMBLY_HEADER_FORMAT {
	container = "docker://metagenlab/diag-pipeline-python-r:1.1"

	tag "$meta.id"

	publishDir "${params.outdir}/${meta.id}/final_assembly", mode: 'copy'

	input:
	tuple val(meta), file(assembly)

	output:
	tuple val(meta), path("./${meta.id}.fasta"), emit: formatted_assembly
	path("./${meta.id}.fasta"), emit: ref_assembly
	shell:
	'''
	gzip -d -c !{assembly} > contigs.fasta
	awk '/^>/{print ">contig_" ++i; next}{print}' < contigs.fasta > !{meta.id}.fasta
	'''
}
