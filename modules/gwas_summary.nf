process GWAS_SUMMARY {
    tag "${label}"
    label 'process_low'

    input:
    tuple val(label), path(results)

    output:
    path "${label}_summary.tsv",   emit: summary
    path "${label}_top_hits.tsv",  emit: top_hits
    path "${label}_*.png",         emit: plots

    script:
    """
    gwas_summary.R ${label} ${results} ${params.p_threshold} ${params.n_top_hits}
    """
}
