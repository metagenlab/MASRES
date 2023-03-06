process REPORT_IDENTIFICATION {
	// docker://metagenlab/diag-pipeline-python-r:1.1
	container =  "metagenlab/python-docutils:latest"
	
	input:
	path(mlst_file)
	path(gtdbtk_summary)
	path(skani)
	path(mash_files)
	path(checkm_result)
    path(centrifuge_files)
	

	output:
	path("report_identification.rst"), emit: rst
	path("report_identification.html"), emit: html
	path("report_identification.tsv"), emit: tsv
	path("centrifuge/*"), emit: centrifuge

	script:
	"""
	identification_report.py -m ${mash_files} -c ${centrifuge_files} -g ${gtdbtk_summary} -k ${checkm_result} -o report_identification -s ${skani}
	"""
}
