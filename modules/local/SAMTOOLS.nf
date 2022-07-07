//managing sam files with SAMTOOLS

process SAMTOOLS{
	container = "docker://staphb/samtools"

	tag "$ID"

	input:
	val(mode)
	tuple val(ID), file(input_sam_file)

	output:
	val ID, emit: ID
	tuple val(ID), file(samtools_output), emit: samtools_out

	script:
	samtools_output = "${mode}_samtools.depth"
	"""
	samtools sort -o aln_sorted ${input_sam_file}
	samtools view -h -o aln_sorted.sam aln_sorted
	samtools depth -a aln_sorted.sam > ${samtools_output}
	"""
} 
