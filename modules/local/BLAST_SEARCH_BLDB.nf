process BLDB_SEARCH {
        container= "docker://quay.io/biocontainers/fasta3:36.3.8--h516909a_5"

	tag "$meta.id"

	label 'process_medium'

        input:
        tuple val(meta), file(annotated_proteins)
	path(BLDB_DB)
	
        output:
        tuple val(meta), file(BLDB_out), emit: bldb_tsv

        script:
	BLDB_out= "BLDB_${meta.id}_raw.tsv"
        """
        ssearch36 -T ${task.cpus} -m 8 -b 10 -d 0 -E 1e-5 ${annotated_proteins} ${BLDB_DB} -z 1 > ${BLDB_out}
        """
}

