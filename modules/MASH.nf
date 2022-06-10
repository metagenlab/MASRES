//Assembly QC with Mash

process MASH {
	container = "docker://quay.io/biocontainers/mash:2.2.1--h3d38be6_0"

	tag "$ID"

	publishDir "${params.outdir}/${ID}/QC/", mode: 'copy'

	input:
	tuple val(ID), file(assembly), file(depth)
	path db

	output:
	val ID, emit: ID
	tuple val(ID), path('./distances.tsv'), emit: mash_results

	script:
	"""
	mash screen -w ${db} ${assembly} > distances.tsv
	"""
}
