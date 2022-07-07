//Illumina QC with FastQC


process FASTQC {
        container = "docker://staphb/fastqc"

        tag "$meta.id"

        publishDir "${params.outdir}/${meta.id}/QC", mode: 'copy'

        input:
        tuple val(meta), file(reads)

        output:
        path "${reads[0].getSimpleName()}_fastqc.html"
	path "${reads[1].getSimpleName()}_fastqc.html"

        script:
        """
        fastqc ${reads[0]} ${reads[1]} --threads $params.threads
        """
}
