// Import modules

include{ NANOPLOT    } from '../modules/NANOPLOT'
include{ FLYE        } from '../modules/FLYE'
include{ MEDAKA      } from '../modules/MEDAKA'
include{ HOMOPOLISH  } from '../modules/HOMOPOLISH'
include{ TRIMMOMATIC } from '../modules/TRIMMOMATIC'
include{ FASTQC      } from '../modules/FASTQC'
include{ BWA         } from '../modules/BWA'
include{ POLYPOLISH  } from '../modules/POLYPOLISH'


workflow HYBRID {
        take:
        input_long
	db_directory
	input_short

        main:

        NANOPLOT(input_long)

        FLYE(input_long)

        input_long
                .join(FLYE.out.assembly)
                .set { ch_for_medaka }

        MEDAKA(ch_for_medaka)

        HOMOPOLISH(MEDAKA.out.medaka_polish, db_directory)

	TRIMMOMATIC(input_short)

        FASTQC(TRIMMOMATIC.out.paired)	

	BWA(HOMOPOLISH.out.homopolished, TRIMMOMATIC.out.paired)

	POLYPOLISH(BWA.out.bwa_aligned)

}
