//Illumina QC with FastQC


process FASTQC {
        container = "docker://staphb/fastqc"

        tag "$ID"

        publishDir "${params.outdir}/${ID}/QC", mode: 'copy'

        input:
        tuple val(ID), file(read1), file(read2)

        output:
        path "${read1.getSimpleName()}_fastqc.html"
	path "${read2.getSimpleName()}_fastqc.html"

        script:
        """
        fastqc ${read1} ${read2} --threads $params.threads
        """
}
