//Formatting rgi tsv file

process RGI_FORMAT {
	container = "docker://metagenlab/diag-pipeline-python-r:1.1"

	tag "$meta.id"
	
	publishDir "${params.outdir}/${meta.id}/resistance", mode: 'copy'

	input:
	tuple val(meta), file(assembly), file(depth), file(fna_file), file(gbff_file), file(rgi_file)

	output:
	tuple val(meta), path("./rgi_${meta.id}_formatted.tsv"), emit: rgi_tsv
	tuple val(meta), path("./CDS_depth.tsv"), emit: CDS_depth
	path("./rgi_${meta.id}_formatted.tsv"), emit: RGI_file

	script:
	"""
	calculate_CDS_depth.py ${fna_file} -g ${gbff_file} -d ${depth}
	rgi_format.py ${gbff_file} -r ${rgi_file} -d CDS_depth.tsv -n ${meta.id} -o rgi_${meta.id}_formatted.tsv
	"""
}

	
	
	
