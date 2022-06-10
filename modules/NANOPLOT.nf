//Nanopore QC with Nanoplot


process NANOPLOT {
        container = "docker://staphb/nanoplot"

	tag "$ID"

        publishDir "${params.outdir}/${ID}/QC", mode: 'copy'

        input:
        tuple val(ID), file(reads)

        output:
        path "${ID}_NanoPlot-report.html"

        script:
        """
        NanoPlot -t $params.threads --fastq ${reads} -p ${ID}_
        """
}
