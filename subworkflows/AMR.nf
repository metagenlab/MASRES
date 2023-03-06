
include{ PLATON                 } from '../modules/local/PLATON'
include{ RGI                    } from '../modules/local/RGI'
include{ RGI_FORMAT             } from '../modules/local/RGI_FORMAT'
include{ BLDB_SEARCH            } from '../modules/local/BLAST_SEARCH_BLDB'
include{ FORMAT_RES             } from '../modules/local/FORMAT_RES'
include{ INDIVIDUAL_REPORT      } from '../modules/local/INDIVIDUAL_REPORT'
include{ INDIVIDUAL_DOCX_REPORT } from '../modules/local/INDIVIDUAL_DOCX_REPORT'
include{ MERGE_RESISTANCE       } from '../modules/local/MERGE_RESISTANCE'



workflow AMR {
	take:
	assembly
	faa 
	gbff
	fna
	mlst
	mash_csv_tuple

	main:

	assembly
		.map{ meta, fasta, depth -> tuple(meta, fasta) }
		.set{ ch_assembly }

	assembly
		.map{ meta, fasta, depth -> tuple(meta, depth) }
		.set{ ch_depth }

	RGI(faa)

	BLDB_SEARCH(faa, params.bldb_db)

	PLATON(ch_assembly, params.platon_db)
	
	ch_rgi_format =	ch_assembly.join(ch_depth)
		.join(fna)
		.join(gbff)
		.join(RGI.out.rgi_tsv)
		.map{ meta, assembly, depth, fna, gbff, tsv -> tuple(meta, assembly, depth, fna, gbff, tsv) }
	
	RGI_FORMAT(ch_rgi_format)
	
	ch_res_format = gbff.join(BLDB_SEARCH.out.bldb_tsv)
		.join(RGI_FORMAT.out.CDS_depth)
		.join(RGI_FORMAT.out.rgi_tsv)
		.join(mlst)
		.join(PLATON.out.plasmid_annot)
		.join(mash_csv_tuple)
	
	
	FORMAT_RES(ch_res_format, params.bldb_db)

	INDIVIDUAL_REPORT(FORMAT_RES.out.merged_channel)

	INDIVIDUAL_DOCX_REPORT(INDIVIDUAL_REPORT.out.rst, params.docx_template)

	
	emit:
	RGI           = RGI_FORMAT.out.RGI_file
	BLDB          = FORMAT_RES.out.BLDB_file
	MERGED	      =	FORMAT_RES.out.merged_res
	MLST          = FORMAT_RES.out.mlst
	PLASMIDS      = FORMAT_RES.out.plasmids


}
