
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
if(params.mode == 'long'){
	input_files = Channel
    			.fromPath(params.reads_csv)
    			.splitCsv(header:true)
    			.map{ row-> tuple(row.sampleId, file(row.readsONT)) }
			.view()
    			
	homopolish_db = Channel
			.fromPath(params.homopolish_db)
			.view()
	}

if(params.mode == 'short'){
	input_files = Channel
    			.fromPath(params.reads_csv)
    			.splitCsv(header:true)
    			.map{ row-> tuple(row.sampleId, file(row.read1), file(row.read2))}
			.view()
	}

if(params.mode == 'hybrid'){
        input_long = Channel
                        .fromPath(params.reads_csv)
                        .splitCsv(header:true)
                        .map{ row-> tuple(row.sampleId, file(row.readsONT)) }
                        .view()

	homopolish_db = Channel
                        .fromPath(params.homopolish_db)
                        .view()


	input_short = Channel
                        .fromPath(params.reads_csv)
                        .splitCsv(header:true)
                        .map{ row-> tuple(row.sampleId, file(row.read1), file(row.read2))}
                        .view()

        }

bakta_db_dir = Channel
                        .fromPath(params.bakta_db)
                        .view()
platon_db_dir = Channel
			.fromPath(params.platon_db)
			.view()


include { LONG   } from '../subworkflows/LONG'
include { SHORT  } from '../subworkflows/SHORT'
include { HYBRID } from '../subworkflows/HYBRID'
include { AMR    } from '../subworkflows/AMR'
workflow flexAMR {

	if (params.mode == 'long'){
        LONG (input_files, homopolish_db)
	ASSEMBLY_OUTPUT = LONG.out.assembly_out
	}

	if (params.mode == 'short'){
	SHORT (input_files)
	ASSEMBLY_OUTPUT = SHORT.out.assembly_out
	}

	if (params.mode == 'hybrid'){
	HYBRID(input_long, homopolish_db, input_short)
	ASSEMBLY_OUTPUT = HYBRID.out.assembly_out
	}
	
	AMR(ASSEMBLY_OUTPUT, bakta_db_dir, platon_db_dir)
}


