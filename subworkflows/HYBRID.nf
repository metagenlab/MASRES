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
include{ DEPTH as DEPTH_SHORT   } from '../modules/local/DEPTH'
include{ DEPTH as DEPTH_LONG    } from '../modules/local/DEPTH'

workflow HYBRID {
        take:
        input_long
	homopolish_db
	input_short

        main:
	
	// QC of raw reads with Nanoplot

        NANOPLOT(input_long)

	// Subsampling, read filtering, assembly and Medaka polishing with Dragonflye

        DRAGONFLYE(input_long)

	// Second round homologous polishing with Homopolish
	
	HOMOPOLISH(DRAGONFLYE.out.contigs, homopolish_db)

	// Trimming short reads

	FASTP(input_short, false, false)

	// QC of trimmed reads

	FASTQC(FASTP.out.reads)

	// Alignment of short reads on polished long read assembly

	BWA(HOMOPOLISH.out.homopolished, FASTP.out.reads)

	// Short read polishing of assembly

	POLYPOLISH(BWA.out.bwa_aligned)
	
	// Format fasta headers

	ASSEMBLY_HEADER_FORMAT(POLYPOLISH.out.hybrid_assembly)

	//Mapping short and long reads to assembly to calculate depth at each position

	short_mode = 'sr'
	DEPTH_SHORT(short_mode, ASSEMBLY_HEADER_FORMAT.out.ref_assembly, FASTP.out.reads)
	
	long_mode = 'map-ont'
	DEPTH_LONG(long_mode, ASSEMBLY_HEADER_FORMAT.out.ref_assembly, input_long)

	//Merging depth files to get coverage of both long and short reads for each position

	HYB_COVERAGE(DEPTH_SHORT.out.depth, DEPTH_LONG.out.depth)
	
	
	MERGED_OUTPUT_CHANNEL = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly.join(HYB_COVERAGE.out.hybrid_coverage)
	
	emit:
	assembly      = ASSEMBLY_HEADER_FORMAT.out.formatted_assembly
        depth         = HYB_COVERAGE.out.hybrid_coverage

	
}
