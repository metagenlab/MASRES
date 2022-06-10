//Polishing with short reads using Polypolish

process POLYPOLISH {
	container = "docker://maximebranger/polypolish"

	tag "$ID"

	publishDir "${params.outdir}/${ID}/04_hybrid_assembly/", mode: 'copy'
		

	input:
	tuple val(ID), file(assembly), file(alignments_1), file(alignments_2)

	output:
	tuple val(ID), file(polished_assembly), emit: hybrid_assembly

	script:
	polished_assembly = "${ID}_hybrid_assembly.fasta"
	"""
	polypolish_insert_filter.py --in1 ${alignments_1} --in2 ${alignments_2} --out1 filtered_1.sam --out2 filtered_2.sam
	polypolish ${assembly} filtered_1.sam filtered_2.sam > ${polished_assembly}
	"""
}
