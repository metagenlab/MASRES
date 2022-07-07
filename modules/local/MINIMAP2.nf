//Sequence mapping using minimap2

process MINIMAP2 {
        container = "https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0"

        tag "$meta.id"

        input:
	val(mode)
        file(assembly)
        tuple val(meta), file(reads)

	output:
	tuple val(meta), file(samtools_output), emit: depth

	script:
	samtools_output = "${mode}_samtools.depth"
	if( mode == 'sr')
		"""
		minimap2 -ax ${mode} ${assembly} ${reads[0]} ${reads[1]} > aln.sam
		samtools sort -o aln_sorted.bam aln.sam	
		samtools view -h -o aln_sorted.sam aln_sorted.bam
		samtools depth -a aln_sorted.sam > ${samtools_output}
		"""
	else
		"""
		minimap2 -ax ${mode} ${assembly} ${reads} > aln.sam
		samtools sort -o aln_sorted.bam aln.sam
                samtools view -h -o aln_sorted.sam aln_sorted.bam
                samtools depth -a aln_sorted.sam > ${samtools_output}
		"""
}
