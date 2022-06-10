//Multiple Sequence alignment using minimap2

process MINIMAP2 {
        container = "docker://staphb/minimap2"

        tag "$ID"

        input:
	val(mode)
        tuple val(ID), file(assembly)
        tuple val(ID), file(read1), file(read2)

	output:
	tuple val(ID), file("./aln.sam"), emit: minimap2_alignment

	script:
	if( mode == 'sr')
		"""
		minimap2 -ax ${mode} ${assembly} ${read1} ${read2} > aln.sam
		"""
	else
		"""
		minimap2 -ax ${mode} ${assembly} ${read1} > aln.sam
		"""
}
