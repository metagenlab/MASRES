include{ MLST                   } from '../modules/nf-core/modules/mlst/main'
include{ BAKTA                  } from '../modules/nf-core/modules/bakta/main'
include{ PLATON                 } from '../modules/local/PLATON'
include{ MASH_SCREEN            } from '../modules/nf-core/modules/mash/screen/main'
include{ RGI                    } from '../modules/local/RGI'
include{ RGI_FORMAT             } from '../modules/local/RGI_FORMAT'
include{ BLDB_SEARCH            } from '../modules/local/BLAST_SEARCH_BLDB'
include{ FORMAT_RES             } from '../modules/local/FORMAT_RES'
include{ MASH_FORMAT            } from '../modules/local/MASH_FORMAT'
include { MERGE_RESISTANCE } from '../modules/local/MERGE_RESISTANCE'



workflow AMR {
	take:
	assembly

	main:

	assembly
		.map{ meta, fasta, depth -> tuple(meta, fasta) }
		.set{ ch_assembly }

	assembly
		.map{ meta, fasta, depth -> tuple(meta, depth) }
		.set{ ch_depth }

	MLST( ch_assembly )

	MASH_SCREEN(ch_assembly, params.mash_db)
	
	MASH_FORMAT(MASH_SCREEN.out.screen)

	BAKTA(ch_assembly, params.bakta_db, [], [])
	
	PLATON(ch_assembly, params.platon_db)

	BAKTA.out.faa
		.map{ meta, prot_faa -> tuple(meta, prot_faa) }
                .set{ ch_bakta_prot }

	RGI(BAKTA.out.faa)

	BLDB_SEARCH(BAKTA.out.faa, params.bldb_db)
	
	ch_rgi_format =	ch_assembly.join(ch_depth)
		.join(BAKTA.out.fna)
		.join(BAKTA.out.gbff)
		.join(RGI.out.rgi_tsv)
		.map{ meta, assembly, depth, fna, gbff, tsv -> tuple(meta, assembly, depth, fna, gbff, tsv) }
	
	RGI_FORMAT(ch_rgi_format)

	ch_res_format = BAKTA.out.gbff.join(BLDB_SEARCH.out.bldb_tsv)
		.join(RGI_FORMAT.out.CDS_depth)
		.join(RGI_FORMAT.out.rgi_tsv)
		.join(MLST.out.tsv)
		.join(PLATON.out.plasmid_annot)

	FORMAT_RES(ch_res_format, params.bldb_db)


	emit:
	RGI           = RGI_FORMAT.out.RGI_file
	BLDB          = FORMAT_RES.out.BLDB_file
	MERGED	      =	FORMAT_RES.out.merged_res
	MLST          = FORMAT_RES.out.mlst
	PLASMIDS      = FORMAT_RES.out.plasmids
	MASH          = MASH_FORMAT.out.formatted_mash

}
