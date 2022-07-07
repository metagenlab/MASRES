//Resistance gene predictions with Resistance Gene Identifier

process RGI {
	container= "docker://metagenlab/rgi:5.2.1-3.2.2"

	publishDir "${params.outdir}/${meta.id}/resistance", mode: 'copy'

	input:
	tuple val(meta), file(annotated_proteins)

	output:
	val meta, emit: meta
	tuple val(meta), path("./rgi_${meta.id}.tsv"), emit: rgi_tsv

	script:
	"""
	rgi main -t protein -i ${annotated_proteins} -n ${params.threads} -o rgi_${meta.id}.json
	rgi tab -i rgi_${meta.id}.json
	mv rgi_${meta.id}.txt rgi_${meta.id}.tsv
	"""
}
