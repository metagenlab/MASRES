// Import modules

include{ NANOPLOT    } from '../modules/local/NANOPLOT'
include{ NANOSTAT   } from '../modules/local/NANOSTAT'
include{ DRAGONFLYE } from '../modules/nf-core/modules/dragonflye/main'
include{ HOMOPOLISH  } from '../modules/local/HOMOPOLISH'
include{ FASTP       } from '../modules/nf-core/modules/fastp/main'
include{ FASTQC      } from '../modules/nf-core/modules/fastqc/main'
include{ QUAST       } from '../modules/local/QUAST'
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main'
include{ BWA         } from '../modules/local/BWA'
include{ POLYPOLISH  } from '../modules/local/POLYPOLISH'
include{ HYB_COVERAGE} from '../modules/local/HYB_COV'
include{ ASSEMBLY_HEADER_FORMAT       } from '../modules/local/ASSEMBLY_HEADER_FORMAT'
include{ DEPTH as DEPTH_SHORT   } from '../modules/local/DEPTH'
include{ DEPTH as DEPTH_LONG    } from '../modules/local/DEPTH'
include { CENTRIFUGE_CENTRIFUGE } from '../modules/nf-core/modules/centrifuge/centrifuge/main'
include { CENTRIFUGE_LONG } from '../modules/local/CENTRIFUGE_LONG'
include { QUALIMAP_BAMQC as QUALIMAP_BAMQC_SHORT } from '../modules/nf-core/modules/qualimap/bamqc/main'
include { QUALIMAP_BAMQC as QUALIMAP_BAMQC_LONG } from '../modules/nf-core/modules/qualimap/bamqc/main'
include { BANDAGE_IMAGE   } from '../modules/local/BANDAGE_IMAGE'

workflow HYBRID {
        take:
        input_long
	input_short

        main:

	ch_multiqc_config        =    file("./assets/multiqc_config.yml", checkIfExists: true)
	
	// QC of raw reads with Nanoplot and Nanostat

        NANOPLOT(input_long)

	NANOSTAT(input_long)


	// Subsampling, read filtering, assembly and Medaka polishing with Dragonflye

        DRAGONFLYE(input_long)

	// Second round homologous polishing with Homopolish
	
	HOMOPOLISH(DRAGONFLYE.out.contigs, params.homopolish_db)

	// Trimming short reads

	FASTP(input_short, false, false)

	// QC of trimmed reads

	FASTQC(FASTP.out.reads)

	// Alignment of short reads on polished long read assembly

	ch_bwa = HOMOPOLISH.out.homopolished.join(FASTP.out.reads)

	BWA(ch_bwa)

	// Short read polishing of assembly

	POLYPOLISH(BWA.out.bwa_aligned)
	
	// Format fasta headers

	ASSEMBLY_HEADER_FORMAT(POLYPOLISH.out.hybrid_assembly)

	// Assembly QC quast

	QUAST (ASSEMBLY_HEADER_FORMAT.out.formatted_assembly, [], false, false)

	// Assembly graph image with Bandage

	BANDAGE_IMAGE(DRAGONFLYE.out.gfa)

	// Contamination check on raw reads

	CENTRIFUGE_CENTRIFUGE(FASTP.out.reads, params.centrifuge_db, false, false, false)

	CENTRIFUGE_LONG(input_long, params.centrifuge_db, false, false, false)
	
	//Mapping short and long reads to assembly to calculate depth at each position

	depth_ch_short = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly
                                        .join(FASTP.out.reads)

	depth_ch_long = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly
                                        .join(input_long)

	short_mode = 'sr'
	DEPTH_SHORT(short_mode, depth_ch_short)
	
	long_mode = 'map-ont'
	DEPTH_LONG(long_mode, depth_ch_long)

	//Create qualimap reports for each bam file

	QUALIMAP_BAMQC_SHORT(DEPTH_SHORT.out.bam, [])
	QUALIMAP_BAMQC_LONG(DEPTH_LONG.out.bam, [])

	//Merging depth files to get coverage of both long and short reads for each position

	ch_hybrid = DEPTH_SHORT.out.depth.join(DEPTH_LONG.out.depth)

	HYB_COVERAGE(ch_hybrid)
	
	//MultiQC report

	ch_multiqc_files = Channel.empty()
	ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1][0]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(NANOSTAT.out.txt.collect())
	ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.results.collect())
	ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC_SHORT.out.results.collect{it[1]}.ifEmpty([]))	
	ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC_LONG.out.results.collect{it[1]}.ifEmpty([]))

	MULTIQC(
		ch_multiqc_files.collect(), [ [ch_multiqc_config], [] ])
	emit:
	assembly      = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly.join(HYB_COVERAGE.out.hybrid_coverage)
	assembly_no_name = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly.map{ meta, assembly -> assembly}

	
}
