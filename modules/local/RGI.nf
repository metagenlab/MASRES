//Resistance gene predictions with Resistance Gene Identifier

process RGI {
	container= "docker://metagenlab/rgi:5.2.1-3.2.2"

	tag "$meta.id"

	label 'process_medium'

	input:
	tuple val(meta), file(annotated_proteins)

	output:
	val meta, emit: meta
	tuple val(meta), path("./rgi_${meta.id}.tsv"), emit: rgi_tsv

	script:
	"""
	rgi main -t protein -i ${annotated_proteins} -n ${task.cpus} -o rgi_${meta.id}.json
	rgi tab -i rgi_${meta.id}.json
	mv rgi_${meta.id}.txt rgi_${meta.id}.tsv
	"""
}
