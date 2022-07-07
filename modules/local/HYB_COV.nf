// merge sam files for hybrid coverage with AWK

process HYB_COVERAGE {
	container = "docker://ubuntu"

	tag "$meta.id"

	input:
	tuple val(meta), file(Depth_long)
	tuple val(meta), file(Depth_short)

	output:
	tuple val(meta.id), file("./merged_coverage.depth"), emit: hybrid_coverage

	script:
	"""
	paste ${Depth_long} ${Depth_short} | awk -v OFS='\t' '{print \$1, \$2, \$3+\$6}' > merged_coverage.depth
	"""
}
