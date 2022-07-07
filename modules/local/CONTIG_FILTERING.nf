//Filtering short contigs


process CONTIG_FILTERING {
        container = "docker://metagenlab/diag-pipeline-python-r:1.1"

        tag "$meta.id"

        publishDir "${params.outdir}/${meta.id}/02_filtered_assembly", mode: 'copy'

        input:
        tuple val(meta), file(assembly)

        output:
        tuple val(meta), path("./${meta.id}_final_assembly.fasta"), emit: assembly
	path("./${meta.id}_final_assembly.fasta"), emit: ref_assembly

        shell:
        '''
        gzip -d -c !{assembly} > contigs.fasta
	awk '/^>/{{print (NR==1)?\$0: \"\\n\" \$0;next}} {{printf \"%s\", \$0}}END{{print \"\"}}' contigs.fasta |  awk \'!/^>/ {{ next }} {{ getline seq }} length(seq) >= 500 {{ print \$0 \"\\n\" seq }}\'  > contigs_500bp.fasta
        sed \"s/NODE_\\([0-9]\\+\\)_.*/contig_\\1/\" contigs_500bp.fasta > !{meta.id}_final_assembly.fasta
        '''
}

