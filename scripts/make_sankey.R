#!/usr/bin/env Rscript
## Up/Down/n.s. category table + Sankey (parallel sets) across comparisons (guide section 4.11).
## Usage: make_sankey.R --degs f1.csv,f2.csv --names n1,n2 --labels l1,l2 --alpha 0.01 --out-table t.tsv --out-pdf p.pdf
this_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(this_dir, "tt_analysis_utils.R"))

a <- parse_cli(commandArgs(trailingOnly = TRUE), required = c("degs", "names", "labels", "alpha", "out_table", "out_pdf"))
files <- split_csv(a$degs); nms <- split_csv(a$names); labs <- split_csv(a$labels); alpha <- as.numeric(a$alpha)
if (!(length(files) == length(nms) && length(nms) == length(labs))) stop("--degs, --names, --labels differ in length", call. = FALSE)

degs <- setNames(lapply(files, read.csv, check.names = FALSE, stringsAsFactors = FALSE), nms)
tab <- build_sankey_table(degs, alpha)          # genes common to all comparisons
dir.create(dirname(a$out_table), recursive = TRUE, showWarnings = FALSE)
write.table(tab, a$out_table, sep = "\t", quote = FALSE, row.names = FALSE)
message(nrow(tab), " genes common to all comparisons")

suppressPackageStartupMessages({ library(ggplot2); library(tidyr) })
cols <- c("Down" = "#66A61E", "n.s." = "grey80", "Up" = "#E72A8A")
dir.create(dirname(a$out_pdf), recursive = TRUE, showWarnings = FALSE)

if (length(nms) == 1) {
  ## a parallel set needs >= 2 axes: with one comparison show the category counts instead
  cnt <- as.data.frame(table(cat = factor(tab[[nms]], levels = names(cols))))
  p <- ggplot(cnt, aes(x = cat, y = Freq, fill = cat)) + geom_col(alpha = 0.8) + geom_text(aes(label = Freq), vjust = -0.3) +
    scale_fill_manual(values = cols, name = "Direction", drop = FALSE) + xlab(labs) + ylab("genes") + theme_classic()
  ggsave(a$out_pdf, p, width = 4.5, height = 4)
} else {
  suppressPackageStartupMessages(library(ggforce))
  ## One row per distinct category combination, weighted by the number of genes (same flows as the
  ## guide's one-row-per-gene input with value = 1, but much faster for ~20k genes).
  combo <- aggregate(list(n = tab$gene_id), tab[nms], length)
  combo$id <- seq_len(nrow(combo))
  combo$first <- factor(combo[[nms[1]]], levels = names(cols))
  long <- pivot_longer(combo, cols = all_of(nms), names_to = "Axis", values_to = "cat")
  long$Axis <- factor(long$Axis, levels = nms, labels = labs)
  long$cat <- factor(long$cat, levels = names(cols))
  mk <- function(variant) {
    p <- ggplot(long, aes(x = Axis, id = id, split = cat, value = n)) +
      geom_parallel_sets(aes(fill = first), alpha = 0.5, axis.width = 0.25)
    p <- if (variant == "guide") {
      p + geom_parallel_sets_axes(axis.width = 0.25, aes(fill = after_stat(label))) +
        geom_parallel_sets_labels(aes(label = after_stat(value)))
    } else {   # fallback if the after_stat() variables are not available in this ggforce version
      p + geom_parallel_sets_axes(axis.width = 0.25, fill = "grey50") + geom_parallel_sets_labels(colour = "white")
    }
    p + scale_fill_manual(values = cols, name = "Direction") + guides(y = "none") + theme_classic() +
      theme(axis.line.x.bottom = element_blank(), axis.ticks.x = element_blank(), panel.grid.major = element_blank(),
            axis.title = element_blank())
  }
  p <- tryCatch({ pp <- mk("guide"); invisible(ggplot_build(pp)); pp },
                error = function(e) { message("guide-style axes failed (", conditionMessage(e), "); using plain axes"); mk("plain") })
  ggsave(a$out_pdf, p, width = 4.5, height = 4)
}
