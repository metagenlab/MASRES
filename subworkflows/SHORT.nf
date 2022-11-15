// Import modules

include{ FASTQC        } from '../modules/nf-core/modules/fastqc/main'
include{ FASTP         } from '../modules/nf-core/modules/fastp/main'
include{ DEPTH         } from '../modules/local/DEPTH'
include{ ASSEMBLY_HEADER_FORMAT } from '../modules/local/ASSEMBLY_HEADER_FORMAT'
include{ UNICYCLER        } from '../modules/nf-core/modules/unicycler/main'
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main'
include { QUAST   } from '../modules/local/QUAST'
include { QUALIMAP_BAMQC } from '../modules/nf-core/modules/qualimap/bamqc/main'
include { CENTRIFUGE_CENTRIFUGE } from '../modules/nf-core/modules/centrifuge/centrifuge/main'


workflow SHORT {
	take:
	input_files
	
	main:
	
	ch_multiqc_config        =    file("./assets/multiqc_config.yml", checkIfExists: true)
	
	// Trim Illumina reads with fastp
	
	FASTP(input_files, false, false)

	// QC post trimming with FastQC

	FASTQC(FASTP.out.reads)

	// Creating channel for spades and short read assembly

	ch_unicycler = FASTP.out.reads.map { meta,reads -> [meta, reads, []]}

	UNICYCLER(ch_unicycler)
	
	// Filtering short contigs

	ASSEMBLY_HEADER_FORMAT(UNICYCLER.out.scaffolds)

	// Assembly quality control

	QUAST (ASSEMBLY_HEADER_FORMAT.out.formatted_assembly, [], false, false)

	// Contamination check with centrifuge

	CENTRIFUGE_CENTRIFUGE(FASTP.out.reads, params.centrifuge_db, false, false, false)

	// Calculating depth for every position in assembled genome

	minimap2_mode = 'sr'

	depth_ch = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly
					.join(FASTP.out.reads)

	DEPTH(minimap2_mode, depth_ch)
	
	QUALIMAP_BAMQC(DEPTH.out.bam, [])

	ch_multiqc_files = Channel.empty()
	ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1][0]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.results.collect())
	ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.results.collect{it[1]}.ifEmpty([]))
	
	MULTIQC(
		ch_multiqc_files.collect(), [ [ch_multiqc_config], [] ])

	emit:
	assembly = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly.join(DEPTH.out.depth)
	assembly_no_name = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly.map{ meta, assembly -> assembly}
}
