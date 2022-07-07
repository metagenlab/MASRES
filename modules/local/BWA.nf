//Sequence mapping using bwa

process BWA {
	container = "docker://pegi3s/bwa"

	tag "$meta.id"

	input:
	tuple val(meta), file(assembly)
	tuple val(meta), file(reads)

	output:
	tuple val(meta), file(assembly), path('./alignments_1.sam'), file('./alignments_2.sam'), emit: bwa_aligned

	script:
	"""
	bwa index ${assembly}
	bwa mem -t ${params.threads} -a ${assembly} ${reads[0]} > alignments_1.sam
	bwa mem -t ${params.threads} -a ${assembly} ${reads[1]} > alignments_2.sam
	"""
}
