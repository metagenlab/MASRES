// Import modules

include{ TRIMMOMATIC   } from '../modules/TRIMMOMATIC'
include{ SPADES        } from '../modules/SPADES'
include{ FASTQC        } from '../modules/FASTQC'

workflow SHORT_ASS {
	take:
	input_files
	
	main:
	
	TRIMMOMATIC(input_files)
	
	FASTQC(TRIMMOMATIC.out.paired)

	SPADES(TRIMMOMATIC.out.paired)

	emit:
	ID        = SPADES.out.ID
	contigs   = SPADES.out.contigs
	graph     = SPADES.out.assembly_fastg
		
}
