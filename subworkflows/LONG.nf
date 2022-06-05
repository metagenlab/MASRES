

// Import modules

include{ NANOPLOT   } from '../modules/NANOPLOT'
include{ FLYE       } from '../modules/FLYE'
include{ MEDAKA     } from '../modules/MEDAKA'
include{ HOMOPOLISH } from '../modules/HOMOPOLISH'
include{ MINIMAP2   } from '../modules/MINIMAP2'
include{ SAMTOOLS   } from '../modules/SAMTOOLS'

workflow LONG {
	take:
	input_files
	homopolish_db
	
	main:
	
	NANOPLOT(input_files)

	FLYE(input_files)
	
	input_files
		.join(FLYE.out.assembly)
		.set { ch_for_medaka }
							
	MEDAKA(ch_for_medaka)
	
	HOMOPOLISH(MEDAKA.out.medaka_polish, homopolish_db)
	
	//adding second unused channel
	input_files
                .join(FLYE.out.assembly)
                .set { ch_for_minimap2 }

	minimap2_mode = "map-ont"

	MINIMAP2(minimap2_mode, HOMOPOLISH.out.homopolished, ch_for_minimap2)	
	
	SAMTOOLS(MINIMAP2.out.minimap2_alignment)

	MERGED_OUTPUT_CHANNEL = HOMOPOLISH.out.homopolished.join(SAMTOOLS.out.samtools_out).view()

	emit:
	assembly_out      = MERGED_OUTPUT_CHANNEL
}
