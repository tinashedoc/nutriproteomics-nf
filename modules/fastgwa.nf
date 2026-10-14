process FASTGWA {
    tag "${label}"
    label 'process_high'

    input:
    tuple val(label), path(pheno)
    tuple val(prefix), path(bed), path(bim), path(fam)
    path grm
    path qcovar
    path covar
    val use_covariates

    output:
    tuple val(label), path("${label}.fastGWA"), emit: results
    path "${label}.log",                         emit: log

    script:
    """
    covar_args=()
    if ${use_covariates}; then
        [[ -s ${qcovar} ]] && covar_args+=(--qcovar ${qcovar})
        [[ -s ${covar} ]] && covar_args+=(--covar ${covar})
    fi

    ${params.gcta} \\
        --bfile ${prefix} \\
        --grm-sparse sparse_grm \\
        --fastGWA-mlm \\
        --pheno ${pheno} \\
        "\${covar_args[@]}" \\
        --maf ${params.gwas_maf} \\
        --geno ${params.gwas_geno} \\
        --thread-num ${task.cpus} \\
        --out ${label}
    """
}
