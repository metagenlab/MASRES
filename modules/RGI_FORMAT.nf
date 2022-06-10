//Formatting rgi tsv file

process RGI_FORMAT {
	container = "docker://metagenlab/diag-pipeline-python-r:1.1"
	
	publishDir "${params.outdir}/${ID}/resistance", mode: 'copy'

	input:
	tuple val(ID), file(assembly), file(depth)
	tuple val(ID), file(fna_file)
	tuple val(ID), file(gbff_file)
	tuple val(ID), file(rgi_file)

	output:
	val ID, emit: ID
	tuple val(ID), path("./rgi_${ID}_formatted.tsv")
	
	script:
	"""
	calculate_CDS_depth.py ${fna_file} -g ${gbff_file} -d ${depth}
	rgi_format.py ${gbff_file} -r ${rgi_file} -d CDS_depth.tsv -n ${ID} -o rgi_${ID}_formatted.tsv
	"""
}

	
	
	
