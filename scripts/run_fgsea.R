#!/usr/bin/env Rscript
## Preranked GSEA with fgsea for one gene-set collection x one comparison (replaces guide section 4.10).
## Usage: run_fgsea.R --rnk x.rnk (from make_rnk.R) --gmt sets.gmt --out x.tsv
##   [--min-size 1] [--max-size 1000000] [--seed 42] [--threads 1]
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(this_dir, "tt_analysis_utils.R"))

a <- parse_cli(commandArgs(trailingOnly = TRUE),
               required = c("rnk", "gmt", "out"),
               defaults = list(min_size = "1", max_size = "1000000", seed = "42", threads = "1"))

ranked <- read.delim(a$rnk, header = FALSE, col.names = c("gene_name", "stat"), quote = "", stringsAsFactors = FALSE)
if (nrow(ranked) == 0) stop("empty ranked list", call. = FALSE)
if (anyDuplicated(ranked$gene_name)) stop("duplicate symbols in the ranked list", call. = FALSE)

suppressPackageStartupMessages(library(fgsea))
pathways <- gmtPathways(a$gmt)
stats <- setNames(ranked$stat, ranked$gene_name)

set.seed(as.integer(a$seed))
nth <- as.integer(a$threads)
bp <- if (nth > 1) BiocParallel::MulticoreParam(workers = nth, RNGseed = as.integer(a$seed)) else BiocParallel::SerialParam(RNGseed = as.integer(a$seed))
## gseaParam = 1 is the 'weighted' scoring scheme (Subramanian et al. 2005), as the guide's -scoring_scheme weighted.
## No nperm: fgsea's adaptive multilevel procedure is used (fgsea >= 1.20), not 10000 fixed permutations.
res <- fgsea(pathways = pathways, stats = stats, minSize = as.integer(a$min_size), maxSize = as.integer(a$max_size),
             gseaParam = 1, BPPARAM = bp)
res <- as.data.frame(res)
res$leadingEdge <- vapply(res$leadingEdge, paste, character(1), collapse = ";")
res <- res[order(res$pval), c("pathway", "size", "ES", "NES", "pval", "padj", "leadingEdge")]
dir.create(dirname(a$out), recursive = TRUE, showWarnings = FALSE)
write.table(res, a$out, sep = "\t", quote = FALSE, row.names = FALSE)
message(nrow(res), " gene sets tested, ", sum(res$padj < 0.05, na.rm = TRUE), " with padj < 0.05")
