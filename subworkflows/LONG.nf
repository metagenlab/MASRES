
// Import modules



include{ NANOPLOT   } from '../modules/local/NANOPLOT'
include{ NANOSTAT   } from '../modules/local/NANOSTAT'
include{ DRAGONFLYE } from '../modules/nf-core/modules/dragonflye/main'
include{ HOMOPOLISH } from '../modules/local/HOMOPOLISH'
include{ DEPTH      } from '../modules/local/DEPTH'                                
include{ ASSEMBLY_HEADER_FORMAT} from '../modules/local/ASSEMBLY_HEADER_FORMAT'
include { QUAST   } from '../modules/local/QUAST'
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main'


workflow LONG {
	take:
	input_files

	main:

	ch_multiqc_config        =    file("./assets/multiqc_config.yml", checkIfExists: true)

	// QC of raw reads with Nanoplot and NanoStat
	
	NANOPLOT(input_files)

	NANOSTAT(input_files)

	// Subsampling, read filtering, assembly and Medaka polishing with Dragonflye

	DRAGONFLYE(input_files)

	// Second round homologous polishing with Homopolish
	
	HOMOPOLISH(DRAGONFLYE.out.contigs, params.homopolish_db)

	// Format fasta headers
	
	ASSEMBLY_HEADER_FORMAT(HOMOPOLISH.out.homopolished)
	
	// Assembly QC with QUAST

	QUAST (ASSEMBLY_HEADER_FORMAT.out.formatted_assembly, [], false, false)
	
	// QC report with Multiqc
	
	ch_multiqc_files = Channel.empty()
	ch_multiqc_files = ch_multiqc_files.mix(NANOSTAT.out.txt.collect())
	ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.results.collect())
	
	MULTIQC(
		ch_multiqc_files.collect(), [ [ch_multiqc_config], [] ])	


	// Calculate sequencing depth for each position

	mapping_mode = "map-ont"	


	depth_ch = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly
                                        .join(input_files)

	
	DEPTH(mapping_mode, depth_ch)	

	emit:
	assembly      = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly.join(DEPTH.out.depth)
	assembly_no_name = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly.map{ meta, assembly -> assembly}
}
