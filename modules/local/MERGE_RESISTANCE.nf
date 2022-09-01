process MERGE_RESISTANCE {
	container = "docker://metagenlab/diag-pipeline-python-r:1.1"
	
	publishDir "${params.outdir}/resistance_report", mode: 'copy'

	input:
	path(tsv_files_merged)
	path(mlst_files)

	output:
	tuple path("./merged_table_resistance.tsv"), emit: res_table
	tuple path("./mlst_summary.tsv"), emit: mlst
	script:

	"""
	merge_resistance_tables.py ${tsv_files_merged}
	cat ${mlst_files} > mlst_summary.tsv
	"""
}
