#!/usr/bin/env Rscript
# Simulate a small test dataset: PLINK text genotypes (sibling pairs plus
# unrelated individuals, 3 chromosomes) and a score/covariate table with one
# causal SNP per score. Run via tests/make_test_data.sh.
set.seed(42)
n_fam <- 150; n_unrel <- 200; snps_per_chr <- 700; chrs <- 1:3
m   <- snps_per_chr * length(chrs)
maf <- runif(m, 0.05, 0.5)

founder <- function() matrix(rbinom(2 * m, 1, rep(maf, each = 2)), nrow = 2)  # 2 haplotypes x m
child   <- function(p1, p2) rbind(p1[cbind(sample(1:2, m, TRUE), 1:m)],
                                  p2[cbind(sample(1:2, m, TRUE), 1:m)])

haps <- list(); fid <- character(0)
for (f in seq_len(n_fam)) {
  p1 <- founder(); p2 <- founder()
  haps <- c(haps, list(child(p1, p2), child(p1, p2)))
  fid  <- c(fid, rep(sprintf("FAM%03d", f), 2))
}
for (u in seq_len(n_unrel)) {
  haps <- c(haps, list(founder()))
  fid  <- c(fid, sprintf("UNR%03d", u))
}
n   <- length(haps)
iid <- sprintf("IND%04d", seq_len(n))
geno <- t(sapply(haps, colSums))                                             # n x m allele counts

alleles <- function(h) ifelse(h == 1, "A", "G")
ped <- cbind(fid, iid, 0, 0, sample(1:2, n, TRUE), -9,
             do.call(cbind, lapply(seq_len(m), function(j)
               cbind(alleles(sapply(haps, `[`, 1, j)), alleles(sapply(haps, `[`, 2, j))))))
write.table(ped, "genotypes.ped", quote = FALSE, row.names = FALSE, col.names = FALSE)
map <- data.frame(chr = rep(chrs, each = snps_per_chr), snp = sprintf("rs%05d", seq_len(m)),
                  cm = 0, pos = rep(seq_len(snps_per_chr) * 10000, length(chrs)))
write.table(map, "genotypes.map", quote = FALSE, row.names = FALSE, col.names = FALSE)

causal <- c(plantdriven = 100, animaldriven = 800, vitaminCsugar = 1500, retinolvitB12 = 2000)
age <- round(rnorm(n, 50, 8))
sex <- as.integer(ped[, 5])
d <- data.frame(FID = fid, IID = iid, age = age, sex = sex)
for (s in names(causal)) {
  g <- geno[, causal[[s]]]
  d[[s]] <- 0.6 * (g - mean(g)) / sd(g) + 0.02 * age + 0.3 * sex + rnorm(n)
}
d$plantdriven[c(3, 10)] <- NA   # exercise missing-value handling
d$age[7] <- NA
write.table(d, "scores.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
