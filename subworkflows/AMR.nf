include{ BAKTA     } from '../modules/BAKTA'
include{ PLATON    } from '../modules/PLATON'

workflow AMR {
	take:
	assembly
	bakta_db
	platon_db

	main:

	BAKTA(assembly, bakta_db)
	
	PLATON(assembly, platon_db)

}
