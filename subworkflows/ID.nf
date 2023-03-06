
include{ GTDBTK_CLASSIFY       } from '../modules/local/gtdbtk_classify'
include{ SKANI                 } from '../modules/local/SKANI'

workflow IDENTIFICATION {
	take:
	assembly_combined


	main:

	
    GTDBTK_CLASSIFY (assembly_combined)
    SKANI(assembly_combined, params.gtdbtk_db_skani)


    emit:
    summary     = GTDBTK_CLASSIFY.out.summary
    ani = SKANI.out.tsv
    
    versions    = GTDBTK_CLASSIFY.out.versions

}
