process SPARSE_GRM {
    tag "${prefix}"
    label 'process_high'

    input:
    tuple val(prefix), path(bed), path(bim), path(fam)
    path snplist

    output:
    path "sparse_grm.grm.*", emit: grm

    script:
    def extract = snplist.name != 'NO_FILE' ? "--extract ${snplist}" : ''
    """
    ${params.gcta} \\
        --bfile ${prefix} \\
        ${extract} \\
        --maf ${params.grm_maf} \\
        --make-grm \\
        --sparse-cutoff ${params.sparse_cutoff} \\
        --thread-num ${task.cpus} \\
        --out sparse_grm
    """
}
