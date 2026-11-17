## Stranded signal plots: one self-contained interactive HTML page (plotly inlined, works from file://).
##   stranded_signal: scripts/stranded_signal.py -> STRANDED_SIGNAL
## Locus tracks (+ strand RNA above the axis, - strand RNA mirrored below, one panel per condition, gene models),
## stranded metagene profiles (sense / antisense relative to gene orientation) and an antisense-fraction table.
## Built when `stranded_plots: true` (default) and the stranded bigWigs are made in this run or already exist
## (see stranded_targets() in common.smk). Works in processing-only runs (no DESeq2 tables: only the loci from the
## config / BED, plus the metagene); skipped, not failed, in analysis-only runs without bigWigs.
## Rerun only this rule after editing loci:  snakemake ... --forcerun stranded_signal  (or just run: the changed
## config values are noticed). Script run via $CONDA_PREFIX/bin/python (the launcher's venv bin/ is first on PATH in
## cluster jobs). Needs pyBigWig and matplotlib (envs/ttpipe_plots.yml).

import shlex


def _as_list(v):
    if v is None or v == "":
        return []
    return [str(x) for x in v] if isinstance(v, (list, tuple)) else [str(v)]


STRANDED_LOCI = _as_list(config.get("track_loci", []))
STRANDED_TOP_DEGS = int(config.get("track_top_degs", 10) or 0)
STRANDED_REGIONS_BED = str(config.get("track_regions_bed", "") or "")
STRANDED_METAGENE_BED = str(config.get("metagene_bed", "") or "")
STRANDED_STATIC = _as_list(config.get("track_static_formats", []))
STRANDED_METAGENE_GROUPS = _as_list(config.get("metagene_groups", ["protein_coding", "lncRNA"]))
_bad_static = [f for f in STRANDED_STATIC if f not in ("svg", "pdf", "png")]
if _bad_static:
    raise ValueError(f"track_static_formats may contain 'svg', 'pdf' and 'png'; got {_bad_static}")
# DESeq2 tables are an input only when they will exist: the analysis runs here and top DEGs are wanted
STRANDED_USE_DEGS = STRANDED_TOP_DEGS > 0 and RUN_ANALYSIS


def stranded_inputs(wildcards):
    d = {
        "bigwigs": expand(BIGWIG_STRAND, condition=CONDITIONS, strand=STRANDS),
        "scale": SCALE_FACTORS,
        "gene_info": GENE_INFO,
        "gtf": GTF,
        "counts_forward": COUNTS_FORWARD,
        "counts_reverse": COUNTS_REVERSE,
        "counts_matrix": COUNTS_MATRIX,
        "samples": config["samplesheet"],      # sample -> condition (the samplesheet; no dependence on the DESeq2 sample table)
    }
    if STRANDED_USE_DEGS:
        d["degs"] = expand(DEGS, comparison=COMPARISON_NAMES)
    if STRANDED_REGIONS_BED:
        d["regions_bed"] = STRANDED_REGIONS_BED
    if STRANDED_METAGENE_BED:
        d["metagene_bed"] = STRANDED_METAGENE_BED
    return d


def _quoted(items):
    return " ".join(shlex.quote(i) for i in items)


rule stranded_signal:
    input:
        unpack(stranded_inputs),
    output:
        STRANDED_SIGNAL,
    log:
        OUT + "/logs/stranded_signal.log",
    benchmark:
        OUT + "/benchmarks/stranded_signal.tsv"
    threads: 2
    conda: ENVS + "/ttpipe_plots.yml"
    params:
        bigwigs = " ".join(shlex.quote(f"{c}={s}={BIGWIG_STRAND.format(condition=c, strand=s)}")
                           for c in CONDITIONS for s in STRANDS),
        condition_order = ",".join(CONDITIONS),
        loci = _quoted(STRANDED_LOCI),
        degs = (" --degs " + " ".join(expand(DEGS, comparison=COMPARISON_NAMES))) if STRANDED_USE_DEGS else "",
        top_degs = STRANDED_TOP_DEGS if STRANDED_USE_DEGS else 0,
        regions_bed = shlex.quote(STRANDED_REGIONS_BED),
        metagene_bed = shlex.quote(STRANDED_METAGENE_BED),
        flank = config.get("track_flank", 5000),
        bin_size = config.get("track_bin_size", 100),
        max_loci = config.get("track_max_loci", 40),
        static = _quoted(STRANDED_STATIC),
        groups = _quoted(STRANDED_METAGENE_GROUPS),
        mflank = config.get("metagene_flank", 3000),
        mmin = config.get("metagene_min_length", 1000),
        strand_min = config.get("strand_table_min_reads", 50),
        min_count = config.get("min_count", 10),
        min_samples = config.get("min_samples", 2),
        banner = shlex.quote(str(config.get("qc_report_banner", "") or "")),
    shell:
        """
        $CONDA_PREFIX/bin/python {SCRIPTS}/stranded_signal.py \
            --bigwig {params.bigwigs} --condition-order {params.condition_order} \
            --scale-factors {input.scale} --gene-info {input.gene_info} --gtf {input.gtf} \
            --counts-forward {input.counts_forward} --counts-reverse {input.counts_reverse} \
            --counts-matrix {input.counts_matrix} --samples {input.samples} \
            --loci {params.loci} {params.degs} --top-degs {params.top_degs} \
            --regions-bed {params.regions_bed} --flank {params.flank} --bin-size {params.bin_size} \
            --max-loci {params.max_loci} --static-formats {params.static} \
            --metagene-groups {params.groups} --metagene-flank {params.mflank} \
            --metagene-min-length {params.mmin} --metagene-bed {params.metagene_bed} \
            --strand-table-min-reads {params.strand_min} --min-count {params.min_count} \
            --min-samples {params.min_samples} --threads {threads} \
            --banner {params.banner} --out {output} 2> {log}
        """
