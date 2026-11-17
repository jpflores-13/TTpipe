## Analysis: DESeq2 (all comparisons, one job) -> fgsea (per collection x comparison) -> Sankey.
## Science follows the TT-seq guide sections 4.9-4.11; deliberate changes are listed in the scripts.
## Rules set only `threads`; memory/time come from the profile (set-resources).

## Sankey axis labels: the treatment condition when all comparisons share one reference
## (the guide's "2H", "6H", "24H"), otherwise the full comparison name.
SANKEY_LABELS = [t for t, r in COMPARISONS] if len({r for _, r in COMPARISONS}) == 1 else COMPARISON_NAMES


## With spike-in size factors (spike_in.normalize has 'deseq2') the variant in rules/spike.smk is used instead.
if not SPIKE_DESEQ2:
    rule deseq2:
        input:
            counts = COUNTS_MATRIX,
            samples = SAMPLES_DESEQ2,
        output:
            degs = expand(DEGS, comparison=COMPARISON_NAMES),
            summary = DEG_SUMMARY,
            norm = NORM_COUNTS,
            vst = VST_COUNTS,
            pca = PCA_PDF,
        params:
            degs_pattern = lambda wildcards: DEGS,   # keeps the literal {comparison} (a lambda switches off wildcard expansion); filled in by the script
            design = config.get("design", "~ condition + replicate"),
            reference = REFERENCE,
            condition_order = ",".join(CONDITIONS),
            comparisons = ",".join(f"{t}:{r}" for t, r in COMPARISONS),
            names = ",".join(COMPARISON_NAMES),
            alpha = config.get("alpha", 0.01),
            lfc = config.get("lfc_threshold", 0),
            min_count = config.get("min_count", 10),
            min_samples = config.get("min_samples", 2),
            script = SCRIPTS + "/run_deseq2.R",
        log:
            OUT + "/logs/deseq2.log"
        benchmark:
            OUT + "/benchmarks/deseq2.txt"
        threads: 2
        conda:
            ENVS + "/ttpipe_r.yml"
        shell:
            """
            Rscript {params.script} \
                --counts {input.counts} --samples {input.samples} \
                --design '{params.design}' --reference {params.reference} \
                --condition-order '{params.condition_order}' \
                --comparisons '{params.comparisons}' --names '{params.names}' \
                --alpha {params.alpha} --lfc-threshold {params.lfc} \
                --min-count {params.min_count} --min-samples {params.min_samples} \
                --degs-pattern '{params.degs_pattern}' --summary {output.summary} \
                --norm {output.norm} --vst {output.vst} --pca {output.pca} > {log} 2>&1
            """


def gsea_gmt(wildcards):
    return GSEA_COLLECTIONS[wildcards.collection]


## Ranked list (symbol <TAB> Wald stat). RNK has no {collection} wildcard, so it cannot be a second
## output of the per-collection fgsea rule (Snakemake requires equal wildcards in all outputs): own rule.
rule rnk:
    input:
        degs = DEGS,
        gene_info = GENE_INFO,
    output:
        RNK,
    params:
        script = SCRIPTS + "/make_rnk.R",
    log:
        OUT + "/logs/rnk_{comparison}.log"
    benchmark:
        OUT + "/benchmarks/rnk_{comparison}.txt"
    threads: 1
    conda:
        ENVS + "/ttpipe_r.yml"
    shell:
        "Rscript {params.script} --degs {input.degs} --gene-info {input.gene_info} --out {output} > {log} 2>&1"


rule fgsea:
    input:
        rnk = RNK,
        gmt = gsea_gmt,
    output:
        GSEA_RESULT,
    params:
        min_size = config.get("gsea_min_size", 1),
        max_size = config.get("gsea_max_size", 1000000),
        seed = config.get("gsea_seed", 42),
        script = SCRIPTS + "/run_fgsea.R",
    log:
        OUT + "/logs/fgsea_{collection}_{comparison}.log"
    benchmark:
        OUT + "/benchmarks/fgsea_{collection}_{comparison}.txt"
    threads: 2
    conda:
        ENVS + "/ttpipe_r.yml"
    shell:
        """
        Rscript {params.script} --rnk {input.rnk} --gmt {input.gmt} --out {output} \
            --min-size {params.min_size} --max-size {params.max_size} --seed {params.seed} \
            --threads {threads} > {log} 2>&1
        """


rule sankey:
    input:
        degs = expand(DEGS, comparison=COMPARISON_NAMES),
    output:
        pdf = SANKEY_PDF,
        table = SANKEY_TABLE,
    params:
        degs = ",".join(expand(DEGS, comparison=COMPARISON_NAMES)),
        names = ",".join(COMPARISON_NAMES),
        labels = ",".join(SANKEY_LABELS),
        alpha = config.get("alpha", 0.01),
        script = SCRIPTS + "/make_sankey.R",
    log:
        OUT + "/logs/sankey.log"
    benchmark:
        OUT + "/benchmarks/sankey.txt"
    threads: 1
    conda:
        ENVS + "/ttpipe_r.yml"
    shell:
        """
        Rscript {params.script} --degs '{params.degs}' --names '{params.names}' --labels '{params.labels}' \
            --alpha {params.alpha} --out-table {output.table} --out-pdf {output.pdf} > {log} 2>&1
        """
