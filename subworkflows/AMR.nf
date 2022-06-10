include{ BAKTA                  } from '../modules/BAKTA'
include{ PLATON                 } from '../modules/PLATON'
include{ QUAST                  } from '../modules/QUAST'
include{ MASH                   } from '../modules/MASH'
include{ RGI                    } from '../modules/RGI'
include{ RGI_FORMAT             } from '../modules/RGI_FORMAT'
include{ ASSEMBLY_HEADER_FORMAT } from '../modules/ASSEMBLY_HEADER_FORMAT'

workflow AMR {
	take:
	assembly
	bakta_db
	platon_db
	mash_db

	main:

	
	QUAST(assembly)

	MASH(assembly, mash_db)
	
	BAKTA(assembly, bakta_db)
	
	PLATON(assembly, platon_db)

	RGI(BAKTA.out.faa)
	
	RGI_FORMAT(assembly, BAKTA.out.fna, BAKTA.out.gbff, RGI.out.rgi_tsv)
}
