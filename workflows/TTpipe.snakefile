## TTpipe: TT-seq (TTchem-seq) processing and analysis.
##
## Processing (run_processing): FASTQ -> FastQC [-> fastp] -> STAR -> MAPQ filter -> Picard dedup
##   -> strand-split BAMs -> replicate-merged BAMs -> bigWig tracks; HTSeq gene counts -> count matrix.
## Analysis (run_analysis): DESeq2 per comparison -> fgsea (GSEA preranked) -> Sankey / PCA.
## Results explorer (results_explorer): one self-contained interactive HTML of the DESeq2 / GSEA results.
## Report (qc_report): one self-contained HTML (sample QC, strandedness, replicate correlation, ...).
## The switches are in config/TTpipeConfig.yaml. Rules are in workflows/rules/.

import os

configfile: "config/TTpipeConfig.yaml"

include: "rules/common.smk"

RUN_NAME = "TTpipe " + os.path.splitext(os.path.basename(config["samplesheet"]))[0]

onstart:
    notify("start", f"Started. Samples: {', '.join(SAMPLES)}. Processing: {RUN_PROCESSING}; "
                    f"analysis: {RUN_ANALYSIS}; report: {RUN_REPORT}.", RUN_NAME)

onsuccess:
    print("TTpipe completed successfully!")
    notify("success", f"Finished successfully.\nResults: {OUT}\nQC report: {QC_REPORT}", RUN_NAME)

onerror:
    notify("error", f"A job failed. Per-job logs: {OUT}/<sample>/logs and {OUT}/logs\n"
                    f"Snakemake log: {log}\n\nLast lines:\n{log_tail(log)}", RUN_NAME)

rule all:
    input:
        all_targets(),

## Sample metadata for DESeq2 (and the report), straight from the samplesheet.
## The samplesheet is an input so that editing it (adding samples, changing conditions) regenerates this table;
## without it Snakemake would keep a stale table from an earlier run (seen on a real run: DESeq2 then failed
## with 'sample names differ between the count matrix and the sample table').
rule sample_table:
    input:
        config["samplesheet"],
    output:
        SAMPLES_DESEQ2,
    run:
        import os
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], "w") as fh:
            fh.write("sample\tcondition\treplicate\n")
            for s in SAMPLES:
                fh.write(f"{s}\t{SAMPLE_TABLE.loc[s, 'Condition']}\t{SAMPLE_TABLE.loc[s, 'Replicate']}\n")

include: "rules/processing.smk"
include: "rules/analysis.smk"
include: "rules/explorer.smk"
include: "rules/report.smk"
include: "rules/multiqc.smk"
