// Import modules

include{ NANOPLOT    } from '../modules/NANOPLOT'
include{ FLYE        } from '../modules/FLYE'
include{ MEDAKA      } from '../modules/MEDAKA'
include{ HOMOPOLISH  } from '../modules/HOMOPOLISH'
include{ TRIMMOMATIC } from '../modules/TRIMMOMATIC'
include{ FASTQC      } from '../modules/FASTQC'
include{ BWA         } from '../modules/BWA'
include{ POLYPOLISH  } from '../modules/POLYPOLISH'
include{ HYB_COVERAGE} from '../modules/HYB_COV'
include{ MINIMAP2 as MINIMAP2_SHORT   } from '../modules/MINIMAP2'
include{ MINIMAP2 as MINIMAP2_LONG    } from '../modules/MINIMAP2'
include{ SAMTOOLS as SAMTOOLS_SHORT   } from '../modules/SAMTOOLS'
include{ SAMTOOLS as SAMTOOLS_LONG    } from '../modules/SAMTOOLS'

workflow HYBRID {
        take:
        input_long
	homopolish_db
	input_short

        main:
	
        NANOPLOT(input_long)

        FLYE(input_long)

        input_long
                .join(FLYE.out.assembly)
                .set { ch_for_medaka }

        MEDAKA(ch_for_medaka)

        HOMOPOLISH(MEDAKA.out.medaka_polish, homopolish_db)

	TRIMMOMATIC(input_short)

        FASTQC(TRIMMOMATIC.out.paired)	

	BWA(HOMOPOLISH.out.homopolished, TRIMMOMATIC.out.paired)

	POLYPOLISH(BWA.out.bwa_aligned)
	
	//Mapping short and long reads to assembly
	short_mode = "sr"
	MINIMAP2_SHORT(short_mode, POLYPOLISH.out.hybrid_assembly, TRIMMOMATIC.out.paired)
	

	long_mode = "map-ont"
	

	input_long
                .join(FLYE.out.assembly)
                .set { ch_for_minimap2_long }

	MINIMAP2_LONG(long_mode, POLYPOLISH.out.hybrid_assembly, ch_for_minimap2_long)
	
	//Calculating depth for short and long reads
	sam_mode_short = "short"
	sam_mode_long = "long"
	SAMTOOLS_SHORT(sam_mode_short, MINIMAP2_SHORT.out.minimap2_alignment)
	SAMTOOLS_LONG(sam_mode_long, MINIMAP2_LONG.out.minimap2_alignment)

	//Merging depth files to get coverage of both long and short reads for each position

	HYB_COVERAGE(SAMTOOLS_LONG.out.samtools_out, SAMTOOLS_SHORT.out.samtools_out)

	MERGED_OUTPUT_CHANNEL = POLYPOLISH.out.hybrid_assembly.join(HYB_COVERAGE.out.hybrid_coverage).view()
	
	emit:
	assembly_out     = MERGED_OUTPUT_CHANNEL
	
}
