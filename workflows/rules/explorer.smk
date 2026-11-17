## Results explorer: one self-contained interactive HTML (volcano / MA, gene search, tables, GSEA tab).
##   results_explorer: scripts/results_explorer.py -> RESULTS_EXPLORER
## Runs after the analysis (also in analysis-only mode) and only when `results_explorer: true`.
## The GSEA files are inputs only when gsea_collections is set. Script run via $CONDA_PREFIX/bin/python
## (the launcher's venv bin/ is first on PATH in cluster jobs). Reuses the ttpipe_qc environment.

import shlex

EXPLORER_GSEA = ([f"{c}={cmp}={GSEA_RESULT.format(collection=c, comparison=cmp)}"
                  for c in GSEA_COLLECTIONS for cmp in COMPARISON_NAMES])


rule results_explorer:
    input:
        degs = expand(DEGS, comparison=COMPARISON_NAMES),
        gene_info = GENE_INFO,
        norm = NORM_COUNTS,
        samples = SAMPLES_DESEQ2,
        summary = DEG_SUMMARY,
        gsea = expand(GSEA_RESULT, collection=list(GSEA_COLLECTIONS), comparison=COMPARISON_NAMES) if GSEA_COLLECTIONS else [],
    output:
        RESULTS_EXPLORER,
    log:
        OUT + "/logs/results_explorer.log",
    benchmark:
        OUT + "/benchmarks/results_explorer.tsv"
    threads: 1
    conda: ENVS + "/ttpipe_qc.yml"
    params:
        alpha = config.get("alpha", 0.01),
        lfc = config.get("lfc_threshold", 0),
        condition_order = ",".join(CONDITIONS),
        banner = shlex.quote(str(config.get("qc_report_banner", "") or "")),
        gsea = " ".join(f"--gsea {shlex.quote(s)}" for s in EXPLORER_GSEA),
    shell:
        """
        $CONDA_PREFIX/bin/python {SCRIPTS}/results_explorer.py --degs {input.degs} \
            --gene-info {input.gene_info} --norm {input.norm} --samples {input.samples} \
            --deg-summary {input.summary} {params.gsea} \
            --alpha {params.alpha} --lfc {params.lfc} --condition-order {params.condition_order} \
            --banner {params.banner} --out {output} 2> {log}
        """
