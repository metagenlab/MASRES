//Assembly quality control with quast

process QUAST{
	container="docker://quay.io/biocontainers/quast:5.0.2--py36pl5321hcac48a8_7"

	tag "$ID"

	publishDir "${params.outdir}/${ID}/QC/", mode: 'copy'

	input:
	tuple val(ID), file(assembly), file(depth)

	output:
	val ID, emit: ID
	tuple val(ID), path('./07_QuastQC/results.txt'), emit: quast_results

	script"
	"""
	quast.py ${assembly} -o 07_QuastQC
	"""

	
