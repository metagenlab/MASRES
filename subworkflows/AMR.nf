
include{ BAKTA                  } from '../modules/nf-core/modules/bakta/main'
include{ PLATON                 } from '../modules/local/PLATON'
include{ QUAST                  } from '../modules/nf-core/modules/quast/main'
include{ MASH_SCREEN            } from '../modules/nf-core/modules/mash/screen/main'
include{ RGI                    } from '../modules/local/RGI'
include{ RGI_FORMAT             } from '../modules/local/RGI_FORMAT'


workflow AMR {
	take:
	assembly
	depth
	bakta_db
	platon_db
	mash_db

	main:

	assembly
		.map{ meta, fasta -> fasta }
		.collect()
		.set{ ch_quast }

	QUAST(ch_quast, [], [], false, false)

	MASH_SCREEN(assembly, mash_db)
	
	BAKTA(assembly, bakta_db, [], [])
	
	PLATON(assembly, platon_db)

	RGI(BAKTA.out.faa)
	
	RGI_FORMAT(assembly, depth, BAKTA.out.fna, BAKTA.out.gbff, RGI.out.rgi_tsv)
}
