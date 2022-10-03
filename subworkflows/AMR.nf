include{ MLST                   } from '../modules/nf-core/modules/mlst/main'
include{ BAKTA                  } from '../modules/nf-core/modules/bakta/main'
include{ PLATON                 } from '../modules/local/PLATON'
include{ SOURMASH_SKETCH        } from '../modules/nf-core/modules/sourmash/sketch/main'
include{ SOURMASH_GATHER        } from '../modules/local/SOURMASH_GATHER'
include{ CHECKM_LINEAGEWF       } from '../modules/nf-core/modules/checkm/lineagewf/main' 
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

	main:

	assembly
		.map{ meta, fasta, depth -> tuple(meta, fasta) }
		.set{ ch_assembly }

	assembly
		.map{ meta, fasta, depth -> tuple(meta, depth) }
		.set{ ch_depth }


	MLST( ch_assembly )

	SOURMASH_SKETCH(ch_assembly)
	
	SOURMASH_GATHER(SOURMASH_SKETCH.out.signatures, params.sourmash_db)

	CHECKM_LINEAGEWF(ch_assembly, "fasta")

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
		.join(SOURMASH_GATHER.out.csv_tuple)
	
	
	FORMAT_RES(ch_res_format, params.bldb_db)

	INDIVIDUAL_REPORT(FORMAT_RES.out.merged_channel)

	INDIVIDUAL_DOCX_REPORT(INDIVIDUAL_REPORT.out.rst, params.docx_template)

	
	emit:
	RGI           = RGI_FORMAT.out.RGI_file
	BLDB          = FORMAT_RES.out.BLDB_file
	MERGED	      =	FORMAT_RES.out.merged_res
	MLST          = FORMAT_RES.out.mlst
	PLASMIDS      = FORMAT_RES.out.plasmids
	MASH          = SOURMASH_GATHER.out.csv

}
