

// Import modules

include{ NANOPLOT   } from '../modules/NANOPLOT'
include{ FLYE       } from '../modules/FLYE'
include{ MEDAKA     } from '../modules/MEDAKA'
include{ HOMOPOLISH } from '../modules/HOMOPOLISH'
include{ MINIMAP2   } from '../modules/MINIMAP2'
include{ SAMTOOLS   } from '../modules/SAMTOOLS'
include{ ASSEMBLY_HEADER_FORMAT} from '../modules/ASSEMBLY_HEADER_FORMAT'

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
	
	ASSEMBLY_HEADER_FORMAT(HOMOPOLISH.out.homopolished)
	
	//adding second unused channel since only one set of reads
	input_files
                .join(FLYE.out.assembly)
                .set { ch_for_minimap2 }

	minimap2_mode = "map-ont"

	MINIMAP2(minimap2_mode, ASSEMBLY_HEADER_FORMAT.out.formatted_assembly, ch_for_minimap2)	
	
	SAMTOOLS(minimap2_mode, MINIMAP2.out.minimap2_alignment)

	MERGED_OUTPUT_CHANNEL = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly.join(SAMTOOLS.out.samtools_out)

	emit:
	assembly_out      = MERGED_OUTPUT_CHANNEL
}
