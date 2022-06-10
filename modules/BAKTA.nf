//Annotation with bakta

process BAKTA {
	container = "docker://quay.io/biocontainers/bakta:1.3.3--pyhdfd78af_0"

	tag "$ID"

	publishDir "${params.outdir}/${ID}/", mode: 'copy'

	input:

	tuple val(ID), file(assembly), file(depth)
	path db

	output:
	val ID, emit: ID
	tuple val(ID), file(depth), emit: depth
	tuple val(ID), path("./05_annotation/${ID}.gff3"), emit: gff3
	tuple val(ID), path("./05_annotation/${ID}.faa"), emit: faa
	tuple val(ID), path("./05_annotation/${ID}.txt"), emit: txt
	tuple val(ID), path("./05_annotation/${ID}.gbff"), emit: gbff
	tuple val(ID), path("./05_annotation/${ID}.ffn"), emit: ffn
	tuple val(ID), path("./05_annotation/${ID}.fna"), emit: fna

	script:
	"""
	bakta --db ${db} --output 05_annotation --prefix ${ID} --locus-tag 'bakta' --threads ${params.threads} --complete ${assembly}
	"""
}
