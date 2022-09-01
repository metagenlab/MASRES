//Polishing with short reads using Polypolish

process POLYPOLISH {
	container = "docker://maximebranger/polypolish"

	tag "$meta.id"

	publishDir "${params.outdir}/${meta.id}/04_hybrid_assembly/", mode: 'copy'
		

	input:
	tuple val(meta), file(assembly), file(alignments_1), file(alignments_2)

	output:
	tuple val(meta), file(polished_assembly), emit: hybrid_assembly

	script:
	polished_assembly = "${meta.id}_hybrid_assembly.fasta.gz"
	"""
	polypolish_insert_filter.py --in1 ${alignments_1} --in2 ${alignments_2} --out1 filtered_1.sam --out2 filtered_2.sam
	polypolish ${assembly} filtered_1.sam filtered_2.sam > ${meta.id}_hybrid_assembly.fasta
	gzip ${meta.id}_hybrid_assembly.fasta
	"""
}
