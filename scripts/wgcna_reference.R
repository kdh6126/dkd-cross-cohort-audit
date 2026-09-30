#!/usr/bin/env Rscript
# Canonical WGCNA hub ranking on one training matrix, as the reference implementation.
#
#   Rscript scripts/wgcna_reference.R expr.tsv y.tsv out.tsv
#
# expr.tsv  samples x genes, header = gene ids, already z-scored within cohort (as the
#           Python comparator receives it).  y.tsv  one label per sample (0/1).
#
# The steps are the WGCNA tutorial's, with its defaults, so that nothing here is our choice:
#   soft power      pickSoftThreshold, first power reaching scale-free R^2 >= 0.80
#                   (powerEstimate); if none does, the power with the largest R^2
#   network         unsigned adjacency, unsigned TOM, one block
#   modules         hierarchical clustering on 1-TOM, cutreeDynamic (deepSplit 2,
#                   minModuleSize 30), mergeCloseModules at 0.25
#   trait module    the module whose eigengene correlates most strongly (|r|) with y
#   hub rank        |kME| = |cor(gene, eigengene)| inside the trait module (signedKME)
# Output: gene, module, kME, in_trait_module, plus power / R^2 / trait-module correlation
# repeated on every row so the Python side can log them.

lib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
.libPaths(c(lib, .libPaths()))
suppressMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 2)

args <- commandArgs(trailingOnly = TRUE)
expr <- read.delim(args[1], check.names = FALSE)
y <- scan(args[2], quiet = TRUE)
datExpr <- as.matrix(expr)
storage.mode(datExpr) <- "double"

gsg <- goodSamplesGenes(datExpr, verbose = 0)
keep <- gsg$goodGenes
datExpr <- datExpr[gsg$goodSamples, keep]
y <- y[gsg$goodSamples]

powers <- c(1:10, seq(12, 30, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "unsigned", verbose = 0)
fit <- sft$fitIndices
r2 <- -sign(fit$slope) * fit$SFT.R.sq
power <- sft$powerEstimate
if (is.na(power)) power <- fit$Power[which.max(r2)]
r2_at <- r2[fit$Power == power]

net <- blockwiseModules(datExpr, power = power, networkType = "unsigned", TOMType = "unsigned",
                        minModuleSize = 30, deepSplit = 2, mergeCutHeight = 0.25,
                        maxBlockSize = ncol(datExpr) + 1, numericLabels = TRUE,
                        pamRespectsDendro = FALSE, verbose = 0)
labels <- net$colors
MEs <- net$MEs
MEs <- MEs[, colnames(MEs) != "ME0", drop = FALSE]
if (ncol(MEs) == 0) stop("no modules found")
me_cor <- sapply(MEs, function(m) cor(m, y))
trait_me <- names(me_cor)[which.max(abs(me_cor))]
trait_mod <- as.integer(sub("ME", "", trait_me))
kme <- signedKME(datExpr, MEs, outputColumnName = "kME")
k_trait <- kme[[paste0("kME", trait_mod)]]

out <- data.frame(gene = colnames(datExpr), module = labels, kME = k_trait,
                  in_trait_module = as.integer(labels == trait_mod),
                  power = power, sft_r2 = r2_at, trait_module_cor = me_cor[[trait_me]],
                  trait_module_size = sum(labels == trait_mod), n_modules = ncol(MEs))
write.table(out, args[3], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("power %d (R^2 %.2f), %d modules, trait module %d (%d genes, r=%.2f)\n",
            power, r2_at, ncol(MEs), trait_mod, sum(labels == trait_mod), me_cor[[trait_me]]))
