// Import modules

include{ TRIMMOMATIC   } from '../modules/TRIMMOMATIC'
include{ SPADES        } from '../modules/SPADES'
include{ FASTQC        } from '../modules/FASTQC'
include{ MINIMAP2      } from '../modules/MINIMAP2'
include{ SAMTOOLS      } from '../modules/SAMTOOLS'

workflow SHORT {
	take:
	input_files
	
	main:
	
	TRIMMOMATIC(input_files)
	
	FASTQC(TRIMMOMATIC.out.paired)

	SPADES(TRIMMOMATIC.out.paired)
	
	minimap2_mode = "sr"	

	MINIMAP2(minimap2_mode, SPADES.out.contigs, TRIMMOMATIC.out.paired)
	
	SAMTOOLS(MINIMAP2.out.minimap2_alignment)

	MERGED_OUTPUT_CHANNEL = SPADES.out.contigs.join(SAMTOOLS.out.samtools_out).view()		

	emit:
	assembly_out      = MERGED_OUTPUT_CHANNEL
}
