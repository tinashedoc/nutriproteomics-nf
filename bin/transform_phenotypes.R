#!/usr/bin/env Rscript
# Transform nutrient pattern proteomic scores into GWAS phenotypes.
#
# For each score: regress on the covariates, inverse rank normalise the
# residuals, and write a GCTA phenotype file (FID IID value, no header).
# Also writes GCTA quantitative (qcovar.txt) and categorical (covar.txt)
# covariate files so the covariates can be included again at the GWAS stage.
#
# Usage:
#   transform_phenotypes.R --input data.tsv \
#     --scores plant=plantdriven,animal=animaldriven \
#     [--fid FID] [--iid IID] [--qcovars age] [--covars sex] \
#     [--n_pcs 10] [--pc_prefix PC] [--outdir .]

parse_args <- function(args) {
  opts <- list(fid = "FID", iid = "IID", qcovars = "age", covars = "sex",
               n_pcs = "10", pc_prefix = "PC", outdir = ".")
  if (length(args) %% 2 != 0) stop("Arguments must be --key value pairs")
  for (i in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[i])
    opts[[key]] <- args[i + 1]
  }
  if (is.null(opts$input) || is.null(opts$scores)) {
    stop("--input and --scores are required")
  }
  opts
}

split_list <- function(x) {
  if (is.null(x) || x == "" || x == "none") character(0) else trimws(strsplit(x, ",")[[1]])
}

# Rank-based inverse normal transformation; NAs are kept in place.
inverse_normal <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

write_gcta <- function(df, path) {
  write.table(df, path, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE, na = "NA")
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

score_pairs <- split_list(opts$scores)
SCORES <- setNames(sub("^[^=]+=", "", score_pairs), sub("=.*$", "", score_pairs))
if (any(!grepl("=", score_pairs)) || any(!grepl("^[A-Za-z0-9_]+$", names(SCORES)))) {
  stop("--scores must be label=column pairs; labels may use letters, digits and _")
}

n_pcs   <- as.integer(opts$n_pcs)
pcs     <- if (n_pcs > 0) paste0(opts$pc_prefix, seq_len(n_pcs)) else character(0)
QCOVARS <- c(split_list(opts$qcovars), pcs)
COVARS  <- split_list(opts$covars)
dir.create(opts$outdir, showWarnings = FALSE, recursive = TRUE)

sep <- if (grepl("\\.csv(\\.gz)?$", opts$input, ignore.case = TRUE)) "," else ""
DATASET <- read.table(opts$input, header = TRUE, sep = sep,
                      stringsAsFactors = FALSE, check.names = FALSE)

required <- c(opts$fid, opts$iid, QCOVARS, COVARS, SCORES)
missing  <- setdiff(required, names(DATASET))
if (length(missing) > 0) {
  stop("Input is missing columns: ", paste(missing, collapse = ", "))
}

ids <- data.frame(FID = DATASET[[opts$fid]], IID = DATASET[[opts$iid]])
for (cv in COVARS) DATASET[[cv]] <- factor(DATASET[[cv]])

log_lines <- character(0)
for (label in names(SCORES)) {
  score   <- SCORES[[label]]
  formula <- reformulate(sprintf("`%s`", c(QCOVARS, COVARS)), response = sprintf("`%s`", score))

  # na.exclude pads residuals with NA so they stay aligned with DATASET rows.
  fit   <- lm(formula, data = DATASET, na.action = na.exclude)
  pheno <- inverse_normal(residuals(fit))

  out <- file.path(opts$outdir, paste0("pheno_", label, ".txt"))
  write_gcta(cbind(ids, pheno = pheno), out)
  line <- sprintf("%s\t%s\t%d\t%.4f", label, score, sum(!is.na(pheno)),
                  summary(fit)$r.squared)
  log_lines <- c(log_lines, line)
  message(line)
}

writeLines(c("label\tscore\tn\tcovariate_r2", log_lines),
           file.path(opts$outdir, "transform_summary.tsv"))
if (length(QCOVARS) > 0) {
  write_gcta(cbind(ids, DATASET[, QCOVARS, drop = FALSE]), file.path(opts$outdir, "qcovar.txt"))
}
if (length(COVARS) > 0) {
  write_gcta(cbind(ids, DATASET[, COVARS, drop = FALSE]), file.path(opts$outdir, "covar.txt"))
}
