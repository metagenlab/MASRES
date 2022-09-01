//convert report to docx

process DOCX_REPORT_CREATION {
	container = "docker://metagenlab/pandoc:2.9.2.1"
	
	publishDir "${params.outdir}/resistance_report", mode: 'copy'

	input:
	path(rst_report)
	path(docx_reference)

	output:
	path("./final_report.docx"), emit: docx

	script:
	"""
	pandoc --reference-doc=${docx_reference} ${rst_report} -o final_report.docx --wrap=preserve
	"""
}

