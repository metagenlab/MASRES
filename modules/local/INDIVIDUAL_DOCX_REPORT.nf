//convert individual report to docx

process INDIVIDUAL_DOCX_REPORT {
	container = "docker://metagenlab/pandoc:2.9.2.1"

	tag "$meta.id"
	
	publishDir "${params.outdir}/resistance_report/individual_reports", mode: 'copy'

	input:
	tuple val(meta), file(rst_report)
	path(docx_reference)

	output:
	path("./${meta.id}_report.docx"), emit: docx

	script:
	"""
	####RMV30MAR2026## pandoc --reference-doc=${docx_reference} ${rst_report} -o ${meta.id}_report.docx --wrap=preserve ###
	
	##ADD## Point exactly to the assets folder in the main pipeline directory
    pandoc --reference-doc=${projectDir}/assets/reference.docx ${rst_report} -o ${meta.id}_report.docx --wrap=preserve
    """
}
