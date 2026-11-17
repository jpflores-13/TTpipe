#!/usr/bin/env Rscript
## DESeq2 for all comparisons in one job (guide section 4.9).
## Usage: run_deseq2.R --counts counts.txt --samples samples.txt --design '~ condition + replicate'
##   --reference DMSO --comparisons '2H:DMSO,6H:DMSO' --names '2H_vs_DMSO,6H_vs_DMSO'
##   --degs-pattern 'dir/degs_{comparison}.csv' --summary f --norm f --vst f --pca f
##   [--condition-order DMSO,2H,6H] [--alpha 0.01] [--lfc-threshold 0] [--min-count 10] [--min-samples 2]
## Optional spike-in normalization (EXPERIMENTAL; only passed when spike_in.enabled is true and 'deseq2' is in spike_in.normalize):
##   --size-factors FILE                    'total' method: sizeFactors(dds) are set from FILE (sample, size_factor), matched by sample NAME,
##                                          instead of estimateSizeFactors()
##   --spike-counts FILE --spike-method genes --control-genes FILE
##                                          'genes' method: the DESeq2 object is built on host + spike-in genes and
##                                          estimateSizeFactors(dds, controlGenes = <spike-in genes>) is run; the spike-in rows are then
##                                          dropped before any test. FILE lists the control gene ids (one per line; default: all spike-in rows)
##   --size-factors-out FILE                the size factors used (sample, size_factor, method)
## The size factors are set BEFORE vst(), so the vst output, the PCA and normalized_counts.tsv all use them.
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(this_dir, "tt_analysis_utils.R"))

a <- parse_cli(commandArgs(trailingOnly = TRUE),
               required = c("counts", "samples", "design", "reference", "comparisons", "names",
                            "degs_pattern", "summary", "norm", "vst", "pca"),
               defaults = list(alpha = "0.01", lfc_threshold = "0", min_count = "10",
                               min_samples = "2", condition_order = "", size_factors = "", spike_counts = "",
                               spike_method = "", control_genes = "", size_factors_out = ""))
alpha <- as.numeric(a$alpha); lfc <- as.numeric(a$lfc_threshold)
min_count <- as.numeric(a$min_count); min_samples <- as.numeric(a$min_samples)

suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })

## ---- inputs: columns matched to the sample table BY NAME ----
counts <- read.delim(a$counts, row.names = 1, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
samples <- read.delim(a$samples, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
if (!all(c("sample", "condition", "replicate") %in% colnames(samples)))
  stop("sample table needs columns sample, condition, replicate", call. = FALSE)
counts <- as.matrix(counts[, match_samples(colnames(counts), samples$sample), drop = FALSE])
if (!is.numeric(counts) || anyNA(counts) || any(counts != round(counts)))
  stop("count matrix must contain only non-missing integer counts", call. = FALSE)
storage.mode(counts) <- "integer"

lv <- condition_levels(samples$condition, a$reference, split_csv(a$condition_order))
cmps <- parse_comparisons(a$comparisons, lv)
nms <- split_csv(a$names)
if (length(nms) != length(cmps)) stop("--names and --comparisons differ in length", call. = FALSE)
design <- validate_design(a$design)

coldata <- data.frame(condition = factor(samples$condition, levels = lv),   # reference level first
                      replicate = factor(samples$replicate), row.names = samples$sample)

## ---- optional spike-in size factors (computed before the host-only object is built) ----
spike_sf <- NULL; sf_method <- "median-of-ratios"
if (nzchar(a$size_factors) && nzchar(a$spike_counts))
  stop("give either --size-factors (method total) or --spike-counts (method genes), not both", call. = FALSE)
if (nzchar(a$size_factors)) {
  tab <- read.delim(a$size_factors, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
  if (!all(c("sample", "size_factor") %in% colnames(tab))) stop("size-factor file needs columns sample, size_factor", call. = FALSE)
  spike_sf <- match_size_factors(tab$sample, tab$size_factor, rownames(coldata)); sf_method <- "spike-in total"
} else if (nzchar(a$spike_counts)) {
  if (a$spike_method != "genes") stop("--spike-counts needs --spike-method genes", call. = FALSE)
  sp <- read.delim(a$spike_counts, row.names = 1, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
  sp <- as.matrix(sp[, match_samples(colnames(sp), samples$sample), drop = FALSE]); storage.mode(sp) <- "integer"
  if (anyDuplicated(c(rownames(counts), rownames(sp)))) stop("spike-in gene ids collide with host gene ids", call. = FALSE)
  ctrl <- if (nzchar(a$control_genes)) readLines(a$control_genes) else rownames(sp)
  ctrl <- ctrl[nzchar(ctrl)]
  full <- DESeqDataSetFromMatrix(countData = rbind(counts, sp), colData = coldata, design = ~ 1)
  full <- estimateSizeFactors(full, controlGenes = control_gene_index(rownames(full), ctrl))
  spike_sf <- sizeFactors(full)[rownames(coldata)]; sf_method <- "spike-in genes"
  message("size factors from ", length(ctrl), " spike-in control gene(s)")
}

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = design)

## pre-filter (guide: >= 10 counts in >= 2 samples)
keep <- rowSums(counts(dds) >= min_count) >= min_samples
message(sum(keep), " of ", length(keep), " genes pass the prefilter (>= ", min_count, " counts in >= ", min_samples, " samples)")
if (!any(keep)) stop("no genes pass the prefilter", call. = FALSE)
dds <- dds[keep, ]
if (!is.null(spike_sf)) sizeFactors(dds) <- spike_sf       # prefilter applied to host genes only; factors set before vst()

## ---- VST + PCA (guide: vst(blind = FALSE), plotPCA by replicate and time) ----
## vst() needs >= 1000 genes (its default nsub); the guide is silent on this. Tiny test data
## fall back to varianceStabilizingTransformation(), which has no such limit.
vsd <- if (nrow(dds) >= 1000) vst(dds, blind = FALSE) else varianceStabilizingTransformation(dds, blind = FALSE)
pc <- plotPCA(vsd, intgroup = c("replicate", "condition"), returnData = TRUE)
pv <- round(100 * attr(pc, "percentVar"))
p <- ggplot(pc, aes(PC1, PC2, colour = condition, shape = replicate)) +
  geom_point(size = 3) + xlab(paste0("PC1: ", pv[1], "% variance")) + ylab(paste0("PC2: ", pv[2], "% variance")) +
  coord_fixed() + theme_bw()
dir.create(dirname(a$pca), recursive = TRUE, showWarnings = FALSE)
ggsave(a$pca, p, width = 5, height = 4)

## ---- DESeq2 pipeline (guide order) ----
if (is.null(spike_sf)) dds <- estimateSizeFactors(dds)    # with spike-in factors they are already set (and used by vst above)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

write_tsv_ids <- function(m, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(data.frame(gene_id = rownames(m), m, check.names = FALSE), path, sep = "\t", quote = FALSE, row.names = FALSE)
}
if (nzchar(a$size_factors_out)) {
  dir.create(dirname(a$size_factors_out), recursive = TRUE, showWarnings = FALSE)
  write.table(data.frame(sample = colnames(dds), size_factor = sizeFactors(dds), method = sf_method),
              a$size_factors_out, sep = "\t", quote = FALSE, row.names = FALSE)
}
write_tsv_ids(counts(dds, normalized = TRUE), a$norm)
write_tsv_ids(assay(vsd), a$vst)

## ---- results per comparison ----
## NOTE: with lfc_threshold > 0, results() tests |log2FC| > threshold (altHypothesis "greaterAbs") and the
## `stat` column is the statistic of THAT test, so ranking by `stat` for GSEA no longer means log2FC/SE.
## The guide uses no threshold for TT-seq (lfc_threshold: 0), where stat is the plain Wald statistic.
summ <- list()
for (i in seq_along(cmps)) {
  res <- results(dds, contrast = c("condition", cmps[[i]][1], cmps[[i]][2]), alpha = alpha, lfcThreshold = lfc)
  df <- as.data.frame(res)
  out <- gsub("{comparison}", nms[i], a$degs_pattern, fixed = TRUE)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(gene_id = rownames(df), df, check.names = FALSE), out, row.names = FALSE, quote = FALSE)
  summ[[i]] <- summarize_deg(nms[i], df$padj, df$log2FoldChange, alpha)
  message(nms[i], ": ", summ[[i]]$n_up, " up, ", summ[[i]]$n_down, " down (padj < ", alpha, ") of ", summ[[i]]$n_tested, " tested")
}
write.table(do.call(rbind, summ), a$summary, sep = "\t", quote = FALSE, row.names = FALSE)
