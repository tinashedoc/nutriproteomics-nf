process TRANSFORM_PHENOTYPES {
    tag "${input.name}"
    label 'process_low'

    input:
    path input

    output:
    path "pheno_*.txt",           emit: phenos
    path "qcovar.txt",            emit: qcovar
    path "covar.txt",             emit: covar
    path "transform_summary.tsv", emit: summary

    script:
    """
    transform_phenotypes.R \\
        --input '${input}' \\
        --scores '${params.scores}' \\
        --fid '${params.fid_col}' \\
        --iid '${params.iid_col}' \\
        --qcovars '${params.qcovars}' \\
        --covars '${params.covars}' \\
        --n_pcs ${params.n_pcs} \\
        --pc_prefix '${params.pc_prefix}' \\
        --outdir .

    # Empty files mean "no covariates of this type"; FASTGWA skips them.
    touch qcovar.txt covar.txt
    """
}
