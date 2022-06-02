//Multiple Sequence alignment using minimap2

process MINIMAP2 {
        container = "docker://staphb/minimap2"

        tag "$ID"

        input:
        tuple val(ID), file(assembly)
        tuple val(ID), file(read1), file(read2)

	output:
	tuple val(ID), file(alignment), emit: minimap2_alignment

	script:
	"""
	minimap2 -ax ${params.mapping_mode_depth} ref.fa read1.fa read2.fa > aln.sam
