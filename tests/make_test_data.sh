#!/usr/bin/env bash
# Regenerate tests/data (needs R and PLINK 1.9 on PATH).
set -euo pipefail
cd "$(dirname "$0")"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

(cd "$work" && Rscript "$OLDPWD/simulate_data.R")
plink --file "$work/genotypes" --make-bed --out "$work/genotypes" >/dev/null
plink --bfile "$work/genotypes" --pca 3 --out "$work/pca" >/dev/null

mkdir -p data
cp "$work"/genotypes.{bed,bim,fam} data/
# Add PC1..PC3 to the score table.
Rscript -e '
  d   <- read.table(commandArgs(TRUE)[1], header = TRUE)
  pcs <- read.table(commandArgs(TRUE)[2])[, 1:5]
  names(pcs) <- c("FID", "IID", "PC1", "PC2", "PC3")
  write.table(merge(d, pcs, sort = FALSE), commandArgs(TRUE)[3], quote = FALSE, sep = "\t", row.names = FALSE)
' "$work/scores.tsv" "$work/pca.eigenvec" data/phenotypes.tsv
echo "Test data written to tests/data"
