// Import modules

include{ SPADES        } from '../modules/nf-core/modules/spades/main'
include{ FASTQC        } from '../modules/nf-core/modules/fastqc/main'
include{ FASTP         } from '../modules/nf-core/modules/fastp/main'
include{ DEPTH         } from '../modules/local/DEPTH'
include{ CONTIG_FILTERING } from '../modules/local/CONTIG_FILTERING'


workflow SHORT {
	take:
	input_files
	
	main:

	// Trim Illumina reads with fastp
	
	FASTP(input_files, false, false)

	// QC post trimming with FastQC

	FASTQC(FASTP.out.reads)

	// Creating channel for spades and short read assembly

	ch_spades = FASTP.out.reads.map { meta,reads -> [meta, reads, [], [] ]}

	SPADES(ch_spades, [])
	
	// Filtering short contigs

	CONTIG_FILTERING(SPADES.out.contigs)

	// Calculating depth for every position in assembled genome

	minimap2_mode = 'sr'

	DEPTH(minimap2_mode, CONTIG_FILTERING.out.ref_assembly, FASTP.out.reads)
	
	emit:
	assembly = CONTIG_FILTERING.out.assembly
	depth = DEPTH.out.depth
}
