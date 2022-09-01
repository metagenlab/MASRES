
process FORMAT_RES {
	container = "docker://metagenlab/diag-pipeline-python-r:1.1"

	tag "$meta.id"
	
	publishDir "${params.outdir}/${meta.id}/resistance", mode: 'copy'

	input:
	tuple val(meta), file(gbff_file), file(BLDB_out), file(CDS_depth), file(rgi_file), file(mlst_file), file(plasmid_file)
	path(BLDB_db)

	output:
	tuple val(meta), path("./BLDB_${meta.id}.tsv"), emit: formatted_BLDB
	path("./BLDB_${meta.id}.tsv"), emit: BLDB_file
	path(out_file), emit: merged_res
	path(mlst_file), emit: mlst
	path(plasmid_file), emit: plasmids

	script:
	out_file="${meta.id}.res"
	"""
	extract_BBH.py ${gbff_file} -b ${BLDB_out} -c ${CDS_depth} -n ${meta.id} -i ${params.id_cutoff} -d ${BLDB_db} -o BLDB_${meta.id}.tsv
	merge_resistance_tables.py ${rgi_file} BLDB_${meta.id}.tsv -o ${out_file}
	"""
}
