// merge sam files for hybrid coverage with AWK

process HYB_COVERAGE {
	container = "docker://ubuntu"

	tag "$ID"

	input:
	tuple val(ID), file(Depth_long)
	tuple val(ID), file(Depth_short)

	output:
	val ID, emit: ID
	tuple val(ID), file("./merged_coverage.depth"), emit: hybrid_coverage

	script:
	"""
	paste ${Depth_long} ${Depth_short} | awk -v OFS='\t' '{print \$1, \$2, \$3+\$6}' > merged_coverage.depth
	"""
}
