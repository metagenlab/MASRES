//Generate docx report

process REPORT_GENERATION {
	container = "docker://metagenlab/diag-pipeline-python-r:1.1"
	
	input:
	path(mlst_file)
	path(rgi_bldb_files)
	path(mash_files)
	path(plasmid_files)
	

	output:
	path("./combined_report.rst"), emit: rst

	script:
	"""
	card_mlst_rst_report.py ${mlst_file} -r ${rgi_bldb_files} -m ${mash_files} -p ${plasmid_files}
	"""
}

