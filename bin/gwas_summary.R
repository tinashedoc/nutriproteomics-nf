#!/usr/bin/env Rscript
# Genomic inflation (lambda), QQ and Manhattan plots and top hits for one
# fastGWA result file.
#
# Usage: gwas_summary.R <label> <results.fastGWA> [p_threshold] [n_top]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: gwas_summary.R <label> <results.fastGWA> [p_threshold] [n_top]")
label  <- args[1]
file   <- args[2]
p_sig  <- if (length(args) >= 3) as.numeric(args[3]) else 5e-8
n_top  <- if (length(args) >= 4) as.integer(args[4]) else 100

res <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
res <- res[!is.na(res$P), ]

lambda <- median(qchisq(res$P, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1)
n_sig  <- sum(res$P < p_sig)
writeLines(c("label\tn_snps\tlambda\tn_significant\tmin_p",
             sprintf("%s\t%d\t%.4f\t%d\t%.3g", label, nrow(res), lambda, n_sig, min(res$P))),
           paste0(label, "_summary.tsv"))

top <- head(res[order(res$P), ], n_top)
write.table(top, paste0(label, "_top_hits.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

png(paste0(label, "_qq.png"), width = 1600, height = 1600, res = 300)
obs <- -log10(sort(res$P))
exp <- -log10(ppoints(length(obs)))
plot(exp, obs, pch = 20, cex = 0.4,
     xlab = expression(Expected ~ -log[10](P)), ylab = expression(Observed ~ -log[10](P)),
     main = sprintf("%s (lambda = %.3f)", label, lambda))
abline(0, 1, col = "red")
invisible(dev.off())

res  <- res[order(res$CHR, res$POS), ]
ends <- tapply(res$POS, res$CHR, max)
offs <- setNames(c(0, cumsum(as.numeric(ends))[-length(ends)]), names(ends))
x    <- res$POS + offs[as.character(res$CHR)]
mids <- offs + ends / 2

png(paste0(label, "_manhattan.png"), width = 3000, height = 1400, res = 300)
plot(x, -log10(res$P), pch = 20, cex = 0.4, xaxt = "n",
     col = c("grey30", "steelblue")[res$CHR %% 2 + 1],
     xlab = "Chromosome", ylab = expression(-log[10](P)), main = label)
axis(1, at = mids, labels = names(ends), cex.axis = 0.6)
abline(h = -log10(p_sig), col = "red", lty = 2)
invisible(dev.off())
