process NANOSTAT {
        container = "docker://quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0"

	tag "$meta.id"

        publishDir "${params.outdir}/${meta.id}/QC", mode: 'copy'

        input:
        tuple val(meta), file(reads)

        output:
        path "${meta.id}_nanostat.txt", emit: txt

        script:
        """
        NanoStat --fastq ${reads} -n "${meta.id}_nanostat.txt"
        """
}
