
// Import modules



include{ NANOPLOT   } from '../modules/local/NANOPLOT'
include{ DRAGONFLYE } from '../modules/nf-core/modules/dragonflye/main'
include{ HOMOPOLISH } from '../modules/local/HOMOPOLISH'
include{ DEPTH      } from '../modules/local/DEPTH'                                
include{ ASSEMBLY_HEADER_FORMAT} from '../modules/local/ASSEMBLY_HEADER_FORMAT'

workflow LONG {
	take:
	input_files
	homopolish_db
	
	main:

	// QC of raw reads with Nanoplot
	
	NANOPLOT(input_files)

	// Subsampling, read filtering, assembly and Medaka polishing with Dragonflye

	DRAGONFLYE(input_files)

	// Second round homologous polishing with Homopolish
	
	HOMOPOLISH(DRAGONFLYE.out.contigs, homopolish_db)

	// Format fasta headers
	
	ASSEMBLY_HEADER_FORMAT(HOMOPOLISH.out.homopolished)
	
	// Calculate sequencing depth for each position

	mapping_mode = "map-ont"	
	
	DEPTH(mapping_mode, ASSEMBLY_HEADER_FORMAT.out.ref_assembly, input_files)	

	emit:
	assembly      = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly
	depth         = DEPTH.out.depth
}
