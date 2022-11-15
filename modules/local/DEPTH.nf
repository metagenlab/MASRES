// Sequence mapping and depth using minimap2 and samtools

process DEPTH {
        container = "https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0"

        tag "$meta.id"

	label 'process_low'

	publishDir "${params.outdir}/${meta.id}/mapping/${mode}", pattern: "*bam", mode: 'copy'

        input:
	val(mode)
        tuple val(meta), file(assembly), file(reads)

	output:
	tuple val(meta), file(samtools_output), emit: depth
	tuple val(meta), file(bam_file), emit: bam

	script:
	samtools_output = "${meta.id}-${mode}-samtools.depth"
	bam_file = "${meta.id}.bam"
	if( mode == 'sr')
		"""
		minimap2 -ax ${mode} -t ${task.cpus} ${assembly} ${reads[0]} ${reads[1]} > aln.sam
		samtools sort -o ${bam_file} aln.sam	
		samtools view -h -o aln_sorted.sam ${bam_file}
		samtools depth -a aln_sorted.sam > ${samtools_output}
		"""
	else
		"""
		minimap2 -ax ${mode} -t ${task.cpus} ${assembly} ${reads} > aln.sam
		samtools sort -o ${bam_file} aln.sam
                samtools view -h -o aln_sorted.sam ${bam_file}
                samtools depth -aa aln_sorted.sam > ${samtools_output}
		"""
}
