
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
    			
	db_directory = Channel
			.fromPath(params.db_dir)
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

	db_directory = Channel
                        .fromPath(params.db_dir)
                        .view()


	input_short = Channel
                        .fromPath(params.reads_csv)
                        .splitCsv(header:true)
                        .map{ row-> tuple(row.sampleId, file(row.read1), file(row.read2))}
                        .view()

        }


include { ONT_ASS   } from '../subworkflows/ONT_assembly'
include { SHORT_ASS } from '../subworkflows/SHORT'
include { HYBRID    } from '../subworkflows/HYBRID'

workflow flexAMR {

	if (params.mode == 'long'){
        ONT_ASS (input_files, db_directory)
	}

	if (params.mode == 'short'){
	SHORT_ASS (input_files)
	}

	if (params.mode == 'hybrid'){
	HYBRID(input_long, db_directory, input_short)
	}
}

