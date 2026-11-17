## Optional MultiQC summary (config `multiqc: true`, processing runs only).
## Written next to TTpipe's own TT-seq report: <outputDir>/qc/TTpipe_multiqc.html.
## MultiQC aggregates the standard tool outputs it has modules for: FastQC, fastp, STAR
## (Log.final.out, ReadsPerGene) and Picard MarkDuplicates. TT-seq-specific checks
## (strandedness, rRNA, replicate correlation, DESeq2 summary) are in TTpipe_qc_report.html.
## MultiQC's HTSeq-count module is NOT relied on here: whether it recognises our per-strand
## count files has not been checked.

def multiqc_inputs(wildcards):
    f = (expand(FASTQC_ZIP, sample=SAMPLES, mate=["1", "2"])
         + expand(STAR_LOG, sample=SAMPLES)
         + expand(STAR_GENECNT, sample=SAMPLES)
         + expand(DEDUP_METRICS, sample=SAMPLES))
    if config.get("trim", False):
        f += expand(FASTP_JSON, sample=SAMPLES)
    return f


rule multiqc:
    input:
        multiqc_inputs,
    output:
        html = MULTIQC_REPORT,
    log:
        OUT + "/logs/multiqc.log"
    benchmark:
        OUT + "/benchmarks/multiqc.tsv"
    threads: 1
    conda: ENVS + "/ttpipe_multiqc.yml"
    params:
        outdir = OUT + "/qc",
        name = "TTpipe_multiqc",
        listing = OUT + "/qc/multiqc_files.txt",
        extra = config.get("multiqc_extra", ""),
    shell:
        """
        mkdir -p {params.outdir}
        printf '%s\\n' {input} > {params.listing}
        multiqc --file-list {params.listing} --outdir {params.outdir} --filename {params.name} \
            --title "TTpipe" --force {params.extra} > {log} 2>&1
        """
