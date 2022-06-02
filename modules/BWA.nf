//Multiple Sequence alignment using bwa

process BWA {
	container = "docker://pegi3s/bwa"

	tag "$ID"

	input:
	tuple val(ID), file(assembly)
	tuple val(ID), file(read1), file(read2)

	output:
	tuple val(ID), file(assembly), path('./alignments_1.sam'), file('./alignments_2.sam'), emit: bwa_aligned

	script:
	"""
	bwa index ${assembly}
	bwa mem -t ${params.threads} -a ${assembly} ${read1} > alignments_1.sam
	bwa mem -t ${params.threads} -a ${assembly} ${read2} > alignments_2.sam
	"""
}
