#!/usr/bin/env Rscript
## Ranked list for preranked GSEA: gene symbol <TAB> DESeq2 Wald statistic, NA removed, sorted descending (guide 4.9).
## Deliberate change vs the guide: ranked by symbol (gene_name), not gene_id, because MSigDB sets use symbols;
## genes without a symbol are dropped, duplicate symbols keep the gene with the largest |stat|.
## NOTE: with lfc_threshold > 0 the `stat` column tests |log2FC| > threshold, so it is no longer log2FC/SE.
## Usage: make_rnk.R --degs degs.csv --gene-info gene_info.tsv --out x.rnk
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(this_dir, "tt_analysis_utils.R"))

a <- parse_cli(commandArgs(trailingOnly = TRUE), required = c("degs", "gene_info", "out"))
degs <- read.csv(a$degs, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("gene_id", "stat") %in% colnames(degs))) stop("DEG table needs gene_id and stat columns", call. = FALSE)
ginfo <- read.delim(a$gene_info, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
if (!all(c("gene_id", "gene_name") %in% colnames(ginfo))) stop("gene info needs gene_id and gene_name columns", call. = FALSE)
ranked <- rank_by_symbol(degs, ginfo)
if (nrow(ranked) == 0) stop("no ranked genes with a symbol and a test statistic", call. = FALSE)
dir.create(dirname(a$out), recursive = TRUE, showWarnings = FALSE)
write.table(ranked, a$out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
message(nrow(ranked), " ranked genes")
