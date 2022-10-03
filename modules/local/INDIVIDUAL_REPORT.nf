//Generate docx report for one sample

process INDIVIDUAL_REPORT {

        container = "docker://metagenlab/diag-pipeline-python-r:1.1"

	tag "$meta.id"

        input:
        tuple val(meta), file(mlst_file), file(rgi_bldb_files), file(mash_files), file(plasmid_files)


        output:
        tuple val(meta), path("./${meta.id}_report.rst"), emit: rst

        script:
        """
        card_mlst_rst_report.py ${mlst_file} -r ${rgi_bldb_files} -m ${mash_files} -p ${plasmid_files} -o ${meta.id}_report.rst
	
        """
}

