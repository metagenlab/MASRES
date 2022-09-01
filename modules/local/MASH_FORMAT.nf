//Format MASH output


process MASH_FORMAT {
        container = "docker://metagenlab/diag-pipeline-python-r:1.1"

        tag "$meta.id"

        publishDir "${params.outdir}/${meta.id}/QC", mode: 'copy'

        input:
        tuple val(meta), file(mash_file)

        output:
        path("./${meta.id}.mash"), emit: formatted_mash
	

        shell:
        '''
	cat !{mash_file} | cut -f6 | sed "s/\\[\\.\\.\\.]//" | sed "s/\\[[0-9]\\+ seqs] //" | cut -f2- -d' ' | paste !{mash_file} - | grep -v "ViralProj" | grep -v "phage" | grep -v "virus" | cut -f1,2,4,7 | sort -gr > !{meta.id}.mash
        '''
}

