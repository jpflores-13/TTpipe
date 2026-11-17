## QC summary table and the static HTML report.
##   qc_summary: scripts/qc_summary.py -> SAMPLE_QC (one row per sample)
##   qc_report : scripts/qc_report.py  -> QC_REPORT (one self-contained HTML file)
## Inputs depend on the run mode (RUN_PROCESSING / RUN_ANALYSIS from common.smk) so the report waits for what
## it will read, and works for processing-only, analysis-only (count matrix only) and full runs. Anything
## missing at run time is reported as skipped, never invented.
## Scripts call $CONDA_PREFIX/bin/python explicitly (the launcher's venv bin/ is first on PATH in cluster jobs).

import shlex

QC_CONFIGFILES = " ".join(shlex.quote(str(f)) for f in workflow.configfiles) or "config/TTpipeConfig.yaml"
# Count matrix the QC reads: the pipeline's own when it is produced here, else the user's analysis_counts.
QC_COUNTS = COUNTS_DESEQ2 if RUN_PROCESSING else COUNTS_MATRIX


def qc_summary_inputs(wildcards):
    # the samplesheet is listed so that editing it reruns the summary (it is read via --samplesheet)
    return [config["samplesheet"]] + _qc_summary_data_inputs(wildcards)


def _qc_summary_data_inputs(wildcards):
    if RUN_PROCESSING:
        return (expand(STAR_LOG, sample=SAMPLES) + expand(STAR_GENECNT, sample=SAMPLES)
                + expand(MAPQ_STATS, sample=SAMPLES) + expand(DEDUP_METRICS, sample=SAMPLES)
                + expand(FINAL_STATS, sample=SAMPLES)
                + [COUNTS_FORWARD, COUNTS_REVERSE, COUNTS_DESEQ2, GENE_INFO])
    if RUN_ANALYSIS:
        # analysis-only: only the count matrix (+ gene table, which the analysis needs anyway)
        return [QC_COUNTS, GENE_INFO]
    # report-only run: use whatever already exists on disk, nothing is built
    return [p for p in [QC_COUNTS, GENE_INFO] if os.path.exists(p)]


def qc_report_inputs(wildcards):
    d = {"sample_qc": SAMPLE_QC, "samplesheet": config["samplesheet"]}
    if RUN_ANALYSIS:
        d["deg"] = DEG_SUMMARY
        d["vst"] = VST_COUNTS
    if RUN_PROCESSING:
        d["scale"] = SCALE_FACTORS
        d["bigwigs"] = expand(BIGWIG, condition=CONDITIONS) + expand(BIGWIG_STRAND, condition=CONDITIONS, strand=STRANDS)
        if config.get("sample_tracks", False):
            d["bigwigs"] += expand(SAMPLE_BIGWIG, sample=SAMPLES) + expand(SAMPLE_BIGWIG_STRAND, sample=SAMPLES, strand=STRANDS)
    return d


rule qc_summary:
    input:
        qc_summary_inputs,
    output:
        SAMPLE_QC,
    log:
        OUT + "/logs/qc_summary.log",
    benchmark:
        OUT + "/benchmarks/qc_summary.tsv"
    threads: 1
    conda: ENVS + "/ttpipe_qc.yml"
    params:
        outdir = OUT,
        samplesheet = config["samplesheet"],
        cfgs = QC_CONFIGFILES,
        counts = QC_COUNTS,
    shell:
        """
        $CONDA_PREFIX/bin/python {SCRIPTS}/qc_summary.py --outdir {params.outdir} --samplesheet {params.samplesheet} \
            --config {params.cfgs} --counts {params.counts} --out {output} 2> {log}
        """


rule qc_report:
    input:
        unpack(qc_report_inputs),
    output:
        QC_REPORT,
    log:
        OUT + "/logs/qc_report.log",
    benchmark:
        OUT + "/benchmarks/qc_report.tsv"
    threads: 2
    conda: ENVS + "/ttpipe_qc.yml"
    params:
        outdir = OUT,
        samplesheet = config["samplesheet"],
        cfgs = QC_CONFIGFILES,
        counts = QC_COUNTS,
        run_processing = int(RUN_PROCESSING),
        run_analysis = int(RUN_ANALYSIS),
        banner = shlex.quote(str(config.get("qc_report_banner", "") or "")),
    shell:
        """
        $CONDA_PREFIX/bin/python {SCRIPTS}/qc_report.py --outdir {params.outdir} --samplesheet {params.samplesheet} \
            --config {params.cfgs} --sample-qc {input.sample_qc} --counts {params.counts} \
            --run-processing {params.run_processing} --run-analysis {params.run_analysis} \
            --banner {params.banner} --out {output} 2> {log}
        """
