//Nanopore QC with Nanoplot


process NANOPLOT {
        container = "docker://staphb/nanoplot"

	tag "$meta.id"

        publishDir "${params.outdir}/${meta.id}/QC", mode: 'copy'

        input:
        tuple val(meta), file(reads)

        output:
        path "${meta.id}_NanoPlot-report.html"

        script:
        """
        NanoPlot -t $params.threads --fastq ${reads} -p ${meta.id}_
        """
}
