include{ MLST                   } from '../modules/nf-core/modules/mlst/main'
include{ BAKTA                  } from '../modules/nf-core/modules/bakta/main'
include{ PLATON                 } from '../modules/local/PLATON'
include{ SOURMASH_SKETCH        } from '../modules/nf-core/modules/sourmash/sketch/main'
include{ SOURMASH_GATHER        } from '../modules/local/SOURMASH_GATHER'
include{ CHECKM_LINEAGEWF       } from '../modules/nf-core/modules/checkm/lineagewf/main' 


workflow ANNOTATE {
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
	
	emit:
	MASH          = SOURMASH_GATHER.out.csv
    MASH_TUPLE    = SOURMASH_GATHER.out.csv_tuple
    GBFF          = BAKTA.out.gbff
    FNA           = BAKTA.out.fna
    FAA           = BAKTA.out.faa
    MLST          = MLST.out.tsv
	CHECKM        = CHECKM_LINEAGEWF.out.checkm_tsv
}
