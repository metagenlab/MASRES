
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
include { IDENTIFICATION } from '../subworkflows/ID'
include { AMR    } from '../subworkflows/AMR'
include { ANNOTATE    } from '../subworkflows/ANNOTATE'
include { MERGE_RESISTANCE } from '../modules/local/MERGE_RESISTANCE'
include { REPORT_GENERATION } from '../modules/local/REPORT_GENERATION'
include { DOCX_REPORT_CREATION} from '../modules/local/DOCX_REPORT_CREATION'
include { REPORT_IDENTIFICATION} from '../modules/local/REPORT_IDENTIFICATION'


workflow flexAMR {
	
	INPUT_CHECK(params.reads_csv)


	if (params.mode == 'long'){
        LONG (INPUT_CHECK.out.longreads)
	ANNOTATE(LONG.out.assembly)
	AMR(LONG.out.assembly, ANNOTATE.out.FAA, ANNOTATE.out.GBFF, ANNOTATE.out.FNA, ANNOTATE.out.MLST, ANNOTATE.out.MASH_TUPLE)
	MERGE_RESISTANCE(AMR.out.MERGED.collect(), AMR.out.MLST.collect())
	REPORT_GENERATION(MERGE_RESISTANCE.out.mlst, AMR.out.MERGED.collect(), ANNOTATE.out.MASH.collect(), AMR.out.PLASMIDS.collect())
	DOCX_REPORT_CREATION(REPORT_GENERATION.out.rst, params.docx_template)
	}

	if (params.mode == 'short'){
	SHORT (INPUT_CHECK.out.shortreads)
	ANNOTATE(SHORT.out.assembly)
	AMR(SHORT.out.assembly, ANNOTATE.out.FAA, ANNOTATE.out.GBFF, ANNOTATE.out.FNA, ANNOTATE.out.MLST, ANNOTATE.out.MASH_TUPLE)
	MERGE_RESISTANCE(AMR.out.MERGED.collect(), AMR.out.MLST.collect())
	REPORT_GENERATION(MERGE_RESISTANCE.out.mlst, AMR.out.MERGED.collect(), ANNOTATE.out.MASH.collect(), AMR.out.PLASMIDS.collect())
	DOCX_REPORT_CREATION(REPORT_GENERATION.out.rst, params.docx_template)

    ANNOTATE.out.FNA.map { meta, fna -> [fna] }.collect()
        .set { ch_id }

	//ch_id.view()

	IDENTIFICATION(ch_id)

	ANNOTATE.out.CHECKM.map { meta, checkm -> [checkm] }.collect().set { ch_checkm }
	// ch_checkm.view()
	SHORT.out.centrifuge.map { meta, centri -> [centri] }.collect().set { ch_centri }
	ch_centri.view()

	REPORT_IDENTIFICATION(MERGE_RESISTANCE.out.mlst, // tsv file  
                          IDENTIFICATION.out.summary, 
						  IDENTIFICATION.out.ani,
						  ANNOTATE.out.MASH.collect(), 
						  ch_checkm, 
						  ch_centri)

	}

	if (params.mode == 'hybrid'){
	HYBRID(INPUT_CHECK.out.longreads, INPUT_CHECK.out.shortreads)
	ANNOTATE(HYBRID.out.assembly)
	AMR(HYBRID.out.assembly, ANNOTATE.out.FAA, ANNOTATE.out.GBFF, ANNOTATE.out.FNA, ANNOTATE.out.MLST, ANNOTATE.out.MASH_TUPLE)
    MERGE_RESISTANCE(AMR.out.MERGED.collect(), AMR.out.MLST.collect())
	REPORT_GENERATION(MERGE_RESISTANCE.out.mlst, AMR.out.MERGED.collect(), ANNOTATE.out.MASH.collect(), AMR.out.PLASMIDS.collect())
	DOCX_REPORT_CREATION(REPORT_GENERATION.out.rst, params.docx_template)
	}

}
