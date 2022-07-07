// Import modules

include{ TRIMMOMATIC   } from '../modules/local/TRIMMOMATIC'
include{ SPADES        } from '../modules/nf-core/modules/spades/main'
include{ FASTQC        } from '../modules/nf-core/modules/fastqc/main'
include{ FASTP         } from '../modules/nf-core/modules/fastp/main'
include{ MINIMAP2      } from '../modules/local/MINIMAP2'
include{ CONTIG_FILTERING } from '../modules/local/CONTIG_FILTERING'


workflow SHORT {
	take:
	input_files
	
	main:
	
	FASTP(input_files, false, false)

	FASTQC(FASTP.out.reads)

	ch_spades = FASTP.out.reads.map { meta,reads -> [meta, reads, [], [] ]}
		.view()


	SPADES(ch_spades, [])

	CONTIG_FILTERING(SPADES.out.contigs)
	
	minimap2_mode = 'sr'

	MINIMAP2(minimap2_mode, CONTIG_FILTERING.out.ref_assembly, FASTP.out.reads)
	
	emit:
	assembly = CONTIG_FILTERING.out.assembly
	depth = MINIMAP2.out.depth
}
