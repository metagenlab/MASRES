process GTDBTK_CLASSIFY {
    tag "-"

    conda (params.enable_conda ? "bioconda::gtdbtk=2.0.0" : null)
    container "ecogenomic/gtdbtk:2.1.1"

    input:
    path("bins/*")

    output:
    path "classify/gtdbtk.*.summary.tsv"               , emit: summary
    path "classify/gtdbtk.*.classify.tree.gz"                   , emit: tree
    path "identify/gtdbtk.bac120.markers_summary.tsv"  , emit: markers
    path "align/gtdbtk.*.msa.fasta.gz"                 , emit: msa
    path "align/gtdbtk.*.user_msa.fasta.gz"                     , emit: user_msa
    path "align/gtdbtk.*.filtered.tsv"                 , emit: filtered
    path "gtdbtk.log"                                  , emit: log
    path "gtdbtk.warnings.log"                         , emit: warnings
    path "identify/gtdbtk.failed_genomes.tsv"          , emit: failed
    path "versions.yml"                                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def pplacer_scratch = params.gtdbtk_pplacer_scratch ? "--scratch_dir pplacer_tmp" : ""
    """
    export GTDBTK_DATA_PATH="${params.gtdbtk_db}"
    if [ ${pplacer_scratch} != "" ] ; then
        mkdir pplacer_tmp
    fi

    gtdbtk classify_wf $args \
                    --genome_dir bins \
                    --prefix "gtdbtk" \
                    --out_dir "\${PWD}" \
                    --cpus ${task.cpus} \
                    --pplacer_cpus ${params.gtdbtk_pplacer_cpus} \
                    ${pplacer_scratch} \
                    --min_perc_aa ${params.gtdbtk_min_perc_aa} \
                    --min_af ${params.gtdbtk_min_af} \
                    --extension fna

    (cd classify; gzip gtdbtk.*.classify.tree)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(gtdbtk --version | sed -n 1p | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
