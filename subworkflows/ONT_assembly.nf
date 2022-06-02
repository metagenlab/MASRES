

// Import modules

include{ NANOPLOT   } from '../modules/NANOPLOT'
include{ FLYE       } from '../modules/FLYE'
include{ MEDAKA     } from '../modules/MEDAKA'
include{ HOMOPOLISH } from '../modules/HOMOPOLISH'

workflow ONT_ASS {
	take:
	input_files
	db_directory
	
	main:
	
	NANOPLOT(input_files)

	FLYE(input_files)
	
	input_files
		.join(FLYE.out.assembly)
		.set { ch_for_medaka }
							
	MEDAKA(ch_for_medaka)
	
	HOMOPOLISH(MEDAKA.out.medaka_polish, db_directory)
	
	emit:
	ID                = HOMOPOLISH.out.ID
	Assembly_polished = HOMOPOLISH.out.homopolished
		
}
