## Pure-logic helpers for the TTpipe analysis scripts. Base R only, so they can be unit-tested
## without DESeq2 / fgsea / ggforce (see tests/test_analysis_scripts.py).

## "--key value" command line -> named list. Keys keep dashes converted to underscores.
parse_cli <- function(args, required = character(), defaults = list()) {
  out <- defaults
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (!grepl("^--", a)) stop("unexpected argument: ", a, call. = FALSE)
    if (i == length(args)) stop("missing value for ", a, call. = FALSE)
    out[[gsub("-", "_", sub("^--", "", a))]] <- args[i + 1]
    i <- i + 2
  }
  miss <- setdiff(required, names(out))
  if (length(miss)) stop("missing required argument(s): ", paste0("--", gsub("_", "-", miss), collapse = ", "), call. = FALSE)
  out
}

split_csv <- function(x) {
  v <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  v[nzchar(v)]
}

## The config `design` may only use columns that exist in the sample table, and must use `condition`
## (the contrasts are on it). Returns the formula.
validate_design <- function(design, allowed = c("condition", "replicate")) {
  f <- tryCatch(stats::as.formula(design), error = function(e)
    stop("design '", design, "' is not a valid R formula: ", conditionMessage(e), call. = FALSE))
  if (length(f) != 2) stop("design '", design, "' must be one-sided, e.g. '~ condition + replicate'", call. = FALSE)
  v <- all.vars(f)
  bad <- setdiff(v, allowed)
  if (length(bad)) stop("design '", design, "' uses unknown variable(s): ", paste(bad, collapse = ", "),
                        "; available: ", paste(allowed, collapse = ", "), call. = FALSE)
  if (!"condition" %in% v) stop("design '", design, "' must include 'condition'", call. = FALSE)
  f
}

## Reorder count-matrix columns to the sample table (by name, never by position).
## Returns the column index vector into `count_colnames`.
match_samples <- function(count_colnames, sample_names) {
  if (anyDuplicated(count_colnames)) stop("duplicated sample names in the count matrix header: ",
                                          paste(unique(count_colnames[duplicated(count_colnames)]), collapse = ", "), call. = FALSE)
  only_tab <- setdiff(sample_names, count_colnames)
  only_cnt <- setdiff(count_colnames, sample_names)
  if (length(only_tab) || length(only_cnt))
    stop("sample names differ between the count matrix and the sample table.\n  in table only: ",
         paste(only_tab, collapse = ", "), "\n  in counts only: ", paste(only_cnt, collapse = ", "), call. = FALSE)
  match(sample_names, count_colnames)
}

## Condition levels: reference first, then `order` (if given) or order of appearance.
condition_levels <- function(conditions, reference, order = character()) {
  lv <- if (length(order)) order else unique(conditions)
  if (!reference %in% lv) stop("reference condition '", reference, "' not among conditions: ", paste(lv, collapse = ", "), call. = FALSE)
  if (!all(unique(conditions) %in% lv)) stop("condition_order misses condition(s): ",
                                              paste(setdiff(unique(conditions), lv), collapse = ", "), call. = FALSE)
  c(reference, setdiff(lv, reference))
}

## "treat:ref,treat2:ref2" -> list of c(treat, ref)
parse_comparisons <- function(x, levels) {
  lapply(split_csv(x), function(p) {
    tr <- strsplit(p, ":", fixed = TRUE)[[1]]
    if (length(tr) != 2 || !all(tr %in% levels)) stop("bad comparison '", p, "' (need treatment:reference among ", paste(levels, collapse = ", "), ")", call. = FALSE)
    tr
  })
}

## Up / Down / n.s. (guide section 4.11): NA padj, or padj >= alpha, or log2FC of 0/NA -> n.s.
classify_deg <- function(padj, log2fc, alpha) {
  sig <- !is.na(padj) & padj < alpha
  out <- rep("n.s.", length(padj))
  out[sig & !is.na(log2fc) & log2fc > 0] <- "Up"
  out[sig & !is.na(log2fc) & log2fc < 0] <- "Down"
  out
}

## One row of the DEG summary. n_tested = genes with a non-NA adjusted p-value.
summarize_deg <- function(comparison, padj, log2fc, alpha) {
  cl <- classify_deg(padj, log2fc, alpha)
  data.frame(comparison = comparison, n_up = sum(cl == "Up"), n_down = sum(cl == "Down"),
             n_tested = sum(!is.na(padj)), stringsAsFactors = FALSE)
}

## degs_list: named list (names = comparison names, in axis order) of data.frames with
## gene_id, log2FoldChange, padj. Genes common to all comparisons, in the order of the first.
build_sankey_table <- function(degs_list, alpha) {
  genes <- Reduce(intersect, lapply(degs_list, function(d) d$gene_id))
  tab <- data.frame(gene_id = genes, stringsAsFactors = FALSE)
  for (nm in names(degs_list)) {
    d <- degs_list[[nm]]
    d <- d[match(genes, d$gene_id), ]
    tab[[nm]] <- classify_deg(d$padj, d$log2FoldChange, alpha)
  }
  tab
}

## Ranked list for preranked GSEA: gene symbol + Wald statistic. NA stat and genes without a symbol
## are dropped; duplicate symbols keep the gene with the largest |stat|; sorted descending by stat.
rank_by_symbol <- function(degs, gene_info) {
  sym <- gene_info$gene_name[match(degs$gene_id, gene_info$gene_id)]
  r <- data.frame(gene_name = sym, stat = degs$stat, stringsAsFactors = FALSE)
  r <- r[!is.na(r$stat) & !is.na(r$gene_name) & nzchar(r$gene_name), , drop = FALSE]
  r <- r[order(-abs(r$stat), r$gene_name), , drop = FALSE]     # name as tie-break: deterministic
  r <- r[!duplicated(r$gene_name), , drop = FALSE]
  r <- r[order(-r$stat, r$gene_name), , drop = FALSE]
  rownames(r) <- NULL
  r
}
