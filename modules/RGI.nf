//Resistance gene predictions with Resistance Gene Identifier

process RGI {
	container= "docker://metagenlab/rgi:5.2.1-3.2.2"

	publishDir "${params.outdir}/${ID}/resistance", mode: 'copy'

	input:
	tuple val(ID), file(annotated_proteins)

	output:
	val ID, emit: ID
	tuple val(ID), path("./rgi_${ID}.tsv"), emit: rgi_tsv

	script:
	"""
	rgi main -t protein -i ${annotated_proteins} -n ${params.threads} -o rgi_${ID}.json
	rgi tab -i rgi_${ID}.json
	mv rgi_${ID}.txt rgi_${ID}.tsv
	"""
}
