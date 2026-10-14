# nutriproteomics-nf

A Nextflow pipeline for genome-wide association studies (GWAS) of nutrient pattern
proteomic scores with fastGWA-MLM from [GCTA](https://yanglab.westlake.edu.cn/software/gcta/).

It automates the analysis described in
[tinashedoc/nutriproteomics](https://github.com/tinashedoc/nutriproteomics):

1. **Transform phenotypes.** Each score is regressed on age, sex and genetic principal
   components, and the residuals are inverse rank normalised.
2. **Sparse GRM.** Common SNPs are LD-pruned with PLINK, and a sparse genetic
   relationship matrix (GRM) is built with GCTA.
3. **fastGWA-MLM.** One GWAS runs per phenotype, in parallel. Covariates are included
   again at this stage.
4. **Summary.** The pipeline reports genomic inflation (lambda) and top hits, and draws
   QQ and Manhattan plots for each phenotype.

## Quick start

```bash
# Requires Nextflow (>= 23.10) plus either conda/mamba, or GCTA, PLINK 1.9 and R on PATH.
nextflow run tinashedoc/nutriproteomics-nf -profile test,conda    # simulated test data (or mamba, micromamba)

nextflow run tinashedoc/nutriproteomics-nf -profile conda \
    --input  scores_and_covariates.tsv \
    --bfile  /dataE/AWIGenGWAS/shared/imputed_data_plink/all_imputed_map_qc \
    --outdir results
```

On a SLURM cluster, add `slurm` to the profile list (`-profile conda,slurm`).
Use `-resume` to continue an interrupted run without repeating finished steps.

If GCTA is already installed and is not called `gcta64` on your PATH, point to it with
`--gcta /home/tinashe/GCTA/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1`.

## Input

**Genotypes** (`--bfile`): a PLINK binary fileset prefix (`.bed/.bim/.fam`). Filter imputed
variants on imputation quality (e.g. INFO/R² ≥ 0.3) before converting them to PLINK.

**Scores and covariates** (`--input`): a tab- or whitespace-separated table (`.csv` is
also accepted), with a header and one row per individual:

| Column | Parameter | Default |
|--------|-----------|---------|
| Family and individual IDs, matching the `.fam` file | `--fid_col`, `--iid_col` | `FID`, `IID` |
| Proteomic scores | `--scores` | `plantdriven`, `animaldriven`, `vitaminCsugar`, `retinolvitB12` |
| Quantitative covariates | `--qcovars` | `age` |
| Categorical covariates | `--covars` | `sex` |
| Genetic PCs `<pc_prefix>1..<pc_prefix><n_pcs>` | `--pc_prefix`, `--n_pcs` | `PC`, `10` |

`--scores` takes comma-separated `label=column` pairs. The label names the output files.

```bash
--scores 'plant=plantdriven,animal=animaldriven,vitc=vitaminCsugar,retinol=retinolvitB12'
```

Lists of covariates are comma-separated. Use `none` for no covariates of a type
(e.g. `--covars none`). Missing values are allowed: individuals with missing data for a
score are dropped from that phenotype only.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Score/covariate table |
| `--bfile` | required | PLINK genotype prefix |
| `--scores` | see above | `label=column` pairs |
| `--fid_col` / `--iid_col` | `FID` / `IID` | ID columns |
| `--qcovars` | `age` | Quantitative covariates |
| `--covars` | `sex` | Categorical covariates |
| `--n_pcs` | `10` | Number of genetic PCs |
| `--pc_prefix` | `PC` | PC column prefix |
| `--grm_snps` | none | SNP list for the GRM; overrides pruning |
| `--prune_grm` | `true` | LD-prune SNPs for the GRM with PLINK |
| `--prune_window` / `--prune_step` / `--prune_r2` | `50` / `5` / `0.2` | PLINK `--indep-pairwise` settings |
| `--grm_maf` | `0.01` | MAF filter for GRM SNPs |
| `--sparse_cutoff` | `0.05` | Relatedness cutoff for the sparse GRM |
| `--gwas_covariates` | `true` | Include covariates again in fastGWA |
| `--gwas_maf` | `0.01` | MAF filter for tested SNPs |
| `--gwas_geno` | `0.05` | Maximum SNP missingness |
| `--p_threshold` | `5e-8` | Significance threshold for summaries and plots |
| `--n_top_hits` | `100` | Number of top SNPs reported per phenotype |
| `--gcta` / `--plink` | `gcta64` / `plink` | Tool executables |
| `--outdir` | `results` | Output directory |

CPU, memory and time limits are set per process label (`process_low`, `process_medium`,
`process_high`) in `nextflow.config`. Override them in your own config with `-c my.config`.

## Output

```
results/
├── phenotypes/     pheno_<label>.txt, qcovar.txt, covar.txt, transform_summary.tsv
├── grm/            sparse_grm.grm.{id,sp}, grm_snps.prune.in
├── gwas/           <label>.fastGWA, <label>.log
├── summary/        gwas_summary.tsv (all phenotypes), <label>_top_hits.tsv,
│                   <label>_qq.png, <label>_manhattan.png
└── pipeline_info/  Nextflow report, timeline and trace
```

## Methods notes

- Residuals come from `lm(..., na.action = na.exclude)`, so they stay aligned with
  individuals when data are missing. Inverse normalisation uses
  `qnorm((rank(x) - 0.5) / n)`.
- Covariates are included again at the GWAS stage. Testing inverse-normalised residuals
  without them can lose power and bias effect estimates (Sofer et al., 2019,
  *Genet Epidemiol* 43:263–275).
- If the genetic variance is not significant, fastGWA-MLM switches to linear regression
  for that phenotype; the `.log` file records this.

## Testing

`-profile test` runs the full pipeline on simulated data in `tests/data`: 150 sibling
pairs and 200 unrelated individuals, 2,100 SNPs on 3 chromosomes, and one causal SNP per
score. Each GWAS should report its causal SNP (`rs00100`, `rs00800`, `rs01500`,
`rs02000`) as the only genome-wide significant hit. GitHub Actions runs this on every
push. To regenerate the data, run `tests/make_test_data.sh`.
