//Trimming and filtering with Trimmomatic

process TRIMMOMATIC {
	container = "docker://quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2"
	
	tag "$ID"

	input:
	tuple val(ID), file(read1), file(read2)
	
	output:
	val ID, emit: ID
	tuple val(ID), file(paired1), file(paired2), emit: paired
	tuple val(ID), file(unpaired1), file(unpaired2), emit: unpaired

	script:
	paired1 = "${ID}_R1_paired.fastq"
	paired2 = "${ID}_R2_paired.fastq"
	unpaired1 = "${ID}_R1_unpaired.fastq"
	unpaired2 = "${ID}_R2_unpaired.fastq"
	
	"""
	trimmomatic PE ${read1} ${read2} $paired1 $unpaired1 $paired2 $unpaired2 ILLUMINACLIP:/usr/local/share/trimmomatic-0.39-2/adapters/${params.adapter_file}:${params.adapter_removal_param1}:${params.adapter_removal_param2}:${params.adapter_removal_param3} SLIDINGWINDOW:${params.sliding_window_size}:${params.sliding_window_quality_threshold} LEADING:${params.minqual} TRAILING:${params.minqual} MINLEN:${params.minlen} -threads $params.threads
	"""
}
