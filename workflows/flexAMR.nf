
def helpMessage() {
	log.info"""
	Usage:
	   nextflow run main.nf [Options]
	Options:
	   --reads_csv	Input csv file with fastq paths of raw reads
	   --outdir	Output directory of files (default "results")
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





include { INPUT_CHECK } from '../subworkflows/input_check'
include { LONG   } from '../subworkflows/LONG'
include { SHORT  } from '../subworkflows/SHORT'
include { HYBRID } from '../subworkflows/HYBRID'
include { AMR    } from '../subworkflows/AMR'
include { MERGE_RESISTANCE } from '../modules/local/MERGE_RESISTANCE'
include { REPORT_GENERATION } from '../modules/local/REPORT_GENERATION'
include { DOCX_REPORT_CREATION} from '../modules/local/DOCX_REPORT_CREATION'

workflow flexAMR {
	
	INPUT_CHECK(params.reads_csv)


	if (params.mode == 'long'){
        LONG (INPUT_CHECK.out.longreads)
	AMR(LONG.out.assembly)
	MERGE_RESISTANCE(AMR.out.MERGED.collect(), AMR.out.MLST.collect())
	REPORT_GENERATION(MERGE_RESISTANCE.out.mlst, AMR.out.MERGED.collect(), AMR.out.MASH.collect(), AMR.out.PLASMIDS.collect())
	DOCX_REPORT_CREATION(REPORT_GENERATION.out.rst, params.docx_template)
	}

	if (params.mode == 'short'){
	SHORT (INPUT_CHECK.out.shortreads)
	AMR(SHORT.out.assembly)
	MERGE_RESISTANCE(AMR.out.MERGED.collect(), AMR.out.MLST.collect())
	REPORT_GENERATION(MERGE_RESISTANCE.out.mlst, AMR.out.MERGED.collect(), AMR.out.MASH.collect(), AMR.out.PLASMIDS.collect())
	DOCX_REPORT_CREATION(REPORT_GENERATION.out.rst, params.docx_template)

	}

	if (params.mode == 'hybrid'){
	HYBRID(INPUT_CHECK.out.longreads, INPUT_CHECK.out.shortreads)
	AMR(HYBRID.out.assembly)
        MERGE_RESISTANCE(AMR.out.MERGED.collect(), AMR.out.MLST.collect())
	REPORT_GENERATION(MERGE_RESISTANCE.out.mlst, AMR.out.MERGED.collect(), AMR.out.MASH.collect(), AMR.out.PLASMIDS.collect())
	DOCX_REPORT_CREATION(REPORT_GENERATION.out.rst, params.docx_template)
	}

}


