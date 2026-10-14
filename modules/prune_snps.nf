process PRUNE_SNPS {
    tag "${prefix}"
    label 'process_medium'

    input:
    tuple val(prefix), path(bed), path(bim), path(fam)

    output:
    path "grm_snps.prune.in", emit: snps

    script:
    """
    ${params.plink} \\
        --bfile ${prefix} \\
        --maf ${params.grm_maf} \\
        --geno ${params.gwas_geno} \\
        --indep-pairwise ${params.prune_window} ${params.prune_step} ${params.prune_r2} \\
        --threads ${task.cpus} \\
        --out grm_snps
    """
}
