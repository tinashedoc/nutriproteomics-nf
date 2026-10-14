#!/usr/bin/env nextflow
/*
 * Nutriproteomics GWAS pipeline
 *
 *   1. Regress each proteomic score on covariates and inverse rank normalise the residuals
 *   2. Build a sparse GRM (optionally from LD-pruned common SNPs)
 *   3. Run fastGWA-MLM (GCTA) for each phenotype
 *   4. Summarise each GWAS: lambda, top hits, QQ and Manhattan plots
 */

include { TRANSFORM_PHENOTYPES } from './modules/transform_phenotypes'
include { PRUNE_SNPS           } from './modules/prune_snps'
include { SPARSE_GRM           } from './modules/sparse_grm'
include { FASTGWA              } from './modules/fastgwa'
include { GWAS_SUMMARY         } from './modules/gwas_summary'

// Command-line values such as `--prune_grm false` arrive as strings.
def flag(value) {
    return value.toString().toLowerCase() in ['true', 'yes', '1']
}

workflow {
    if (!params.input) error "Missing --input: table of proteomic scores and covariates"
    if (!params.bfile) error "Missing --bfile: PLINK genotype prefix (without .bed/.bim/.fam)"
    if (!params.scores) error "Missing --scores: label=column pairs, e.g. plant=plantdriven"

    log.info """
        Nutriproteomics GWAS
        ====================
        input     : ${params.input}
        bfile     : ${params.bfile}
        scores    : ${params.scores}
        covariates: ${params.qcovars} + ${params.n_pcs} PCs (quantitative), ${params.covars} (categorical)
        outdir    : ${params.outdir}
        """.stripIndent()

    def bfile = file(params.bfile)
    def genotypes = channel.value(tuple(
        bfile.name,
        file("${params.bfile}.bed", checkIfExists: true),
        file("${params.bfile}.bim", checkIfExists: true),
        file("${params.bfile}.fam", checkIfExists: true)
    ))

    TRANSFORM_PHENOTYPES(file(params.input, checkIfExists: true))

    def phenos = TRANSFORM_PHENOTYPES.out.phenos
        .flatten()
        .map { f -> tuple(f.baseName.replaceFirst(/^pheno_/, ''), f) }

    // SNPs for the GRM: user-supplied list, LD-pruned common SNPs, or all SNPs.
    def grm_snps
    if (params.grm_snps) {
        grm_snps = channel.value(file(params.grm_snps, checkIfExists: true))
    } else if (flag(params.prune_grm)) {
        grm_snps = PRUNE_SNPS(genotypes).snps
    } else {
        grm_snps = channel.value(file("${projectDir}/assets/NO_FILE"))
    }

    SPARSE_GRM(genotypes, grm_snps)

    FASTGWA(
        phenos,
        genotypes,
        SPARSE_GRM.out.grm,
        TRANSFORM_PHENOTYPES.out.qcovar,
        TRANSFORM_PHENOTYPES.out.covar,
        flag(params.gwas_covariates)
    )

    GWAS_SUMMARY(FASTGWA.out.results)

    GWAS_SUMMARY.out.summary
        .collectFile(name: 'gwas_summary.tsv', keepHeader: true, sort: true, storeDir: "${params.outdir}/summary")
}
