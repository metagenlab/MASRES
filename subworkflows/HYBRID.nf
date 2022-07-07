// Import modules

include{ NANOPLOT    } from '../modules/local/NANOPLOT'
include{ DRAGONFLYE } from '../modules/nf-core/modules/dragonflye/main'
include{ HOMOPOLISH  } from '../modules/local/HOMOPOLISH'
include{ FASTP       } from '../modules/nf-core/modules/fastp/main'
include{ FASTQC      } from '../modules/nf-core/modules/fastqc/main'
include{ BWA         } from '../modules/local/BWA'
include{ POLYPOLISH  } from '../modules/local/POLYPOLISH'
include{ HYB_COVERAGE} from '../modules/local/HYB_COV'
include{ ASSEMBLY_HEADER_FORMAT       } from '../modules/local/ASSEMBLY_HEADER_FORMAT'
include{ MINIMAP2 as MINIMAP2_SHORT   } from '../modules/local/MINIMAP2'
include{ MINIMAP2 as MINIMAP2_LONG    } from '../modules/local/MINIMAP2'

workflow HYBRID {
        take:
        input_long
	homopolish_db
	input_short

        main:
	
        NANOPLOT(input_long)

        DRAGONFLYE(input_long)
	
	HOMOPOLISH(DRAGONFLYE.out.contigs, homopolish_db)

	FASTP(input_short, false, false)

	FASTQC(FASTP.out.reads)

	BWA(HOMOPOLISH.out.homopolished, FASTP.out.reads)

	POLYPOLISH(BWA.out.bwa_aligned)
	
	ASSEMBLY_HEADER_FORMAT(POLYPOLISH.out.hybrid_assembly)

	//Mapping short and long reads to assembly

	short_mode = 'sr'
	MINIMAP2_SHORT(short_mode, ASSEMBLY_HEADER_FORMAT.out.ref_assembly, FASTP.out.reads)
	
	long_mode = 'map-ont'
	MINIMAP2_LONG(long_mode, ASSEMBLY_HEADER_FORMAT.out.ref_assembly, input_long)

	//Merging depth files to get coverage of both long and short reads for each position

	HYB_COVERAGE(MINIMAP2_SHORT.out.depth, MINIMAP2_LONG.out.depth)
	
	
	MERGED_OUTPUT_CHANNEL = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly.join(HYB_COVERAGE.out.hybrid_coverage)
	
	emit:
	assembly      = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly
        depth         = HYB_COVERAGE.out.hybrid_coverage

	
}
