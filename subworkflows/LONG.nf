
// Import modules



include{ NANOPLOT   } from '../modules/local/NANOPLOT'
include{ DRAGONFLYE } from '../modules/nf-core/modules/dragonflye/main'
include{ HOMOPOLISH } from '../modules/local/HOMOPOLISH'
include{ MINIMAP2   } from '../modules/local/MINIMAP2'                                
include{ ASSEMBLY_HEADER_FORMAT} from '../modules/local/ASSEMBLY_HEADER_FORMAT'

workflow LONG {
	take:
	input_files
	homopolish_db
	
	main:
	
	NANOPLOT(input_files)

	DRAGONFLYE(input_files)
	
	HOMOPOLISH(DRAGONFLYE.out.contigs, homopolish_db)
	
	ASSEMBLY_HEADER_FORMAT(HOMOPOLISH.out.homopolished)
	
	mapping_mode = "map-ont"	
	
	MINIMAP2(mapping_mode, ASSEMBLY_HEADER_FORMAT.out.ref_assembly, input_files)	

	emit:
	assembly      = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly
	depth         = MINIMAP2.out.depth
}
