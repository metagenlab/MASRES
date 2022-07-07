
def helpMessage() {
	log.info"""
	Usage:
	   nextflow run main.nf [Options]
	Options:
	   --reads_csv	Input csv file with fastq paths of raw reads
	   --outdir	Output directory of files (default "results")
	   --threads    Number of threads to use
	   --mode       Mode to use (short, long or hybrid)
	""".stripIndent()
}

params.help = false

if (params.help) {
    helpMessage()
    exit 0
}



params.outdir = "results"

// Check input

if(!(params.mode in ['long', 'short', 'hybrid'])){
	exit 1, 'Mode selected does not exist, only short, long and hybrid'
	}
if(params.mode == 'long' || params.mode == 'hybrid'){

	homopolish_db = Channel
			.fromPath(params.homopolish_db)
			.view()
	}



bakta_db_dir = Channel
                        .fromPath(params.bakta_db)
                        .view()
platon_db_dir = Channel
			.fromPath(params.platon_db)
			.view()

mash_db_dir = Channel
			.fromPath(params.mash_db)
			.view()

include { INPUT_CHECK } from '../subworkflows/input_check'
include { LONG   } from '../subworkflows/LONG'
include { SHORT  } from '../subworkflows/SHORT'
include { HYBRID } from '../subworkflows/HYBRID'
include { AMR    } from '../subworkflows/AMR'
workflow flexAMR {
	
	INPUT_CHECK(params.reads_csv)


	if (params.mode == 'long'){
        LONG (INPUT_CHECK.out.longreads, homopolish_db)
	ASSEMBLY_OUTPUT = LONG.out.assembly
	DEPTH_OUTPUT = LONG.out.depth
	}

	if (params.mode == 'short'){
	SHORT (INPUT_CHECK.out.shortreads)
	ASSEMBLY_OUTPUT = SHORT.out.assembly
	DEPTH_OUTPUT = SHORT.out.depth
	}

	if (params.mode == 'hybrid'){
	HYBRID(INPUT_CHECK.out.longreads, homopolish_db, INPUT_CHECK.out.shortreads)
	ASSEMBLY_OUTPUT = HYBRID.out.assembly
	DEPTH_OUTPUT = HYBRID.out.depth
	}
	
	AMR(ASSEMBLY_OUTPUT, DEPTH_OUTPUT, bakta_db_dir, platon_db_dir, mash_db_dir)
}


