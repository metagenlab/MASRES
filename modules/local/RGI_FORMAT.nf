//Formatting rgi tsv file

process RGI_FORMAT {
	container = "docker://metagenlab/diag-pipeline-python-r:1.1"
	
	publishDir "${params.outdir}/${meta.id}/resistance", mode: 'copy'

	input:
	tuple val(meta), file(assembly)
	tuple val(meta), file(depth)
	tuple val(meta), file(fna_file)
	tuple val(meta), file(gbff_file)
	tuple val(meta), file(rgi_file)

	output:
	val meta, emit: meta
	tuple val(meta), path("./rgi_${meta.id}_formatted.tsv")
	
	script:
	"""
	calculate_CDS_depth.py ${fna_file} -g ${gbff_file} -d ${depth}
	rgi_format.py ${gbff_file} -r ${rgi_file} -d CDS_depth.tsv -n ${meta.id} -o rgi_${meta.id}_formatted.tsv
	"""
}

	
	
	
