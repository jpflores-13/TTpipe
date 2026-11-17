## Optional spike-in normalization (EXPERIMENTAL). Only included when spike_in.enabled is true (see
## common.smk for the config block and its parse-time checks). With it off, none of these rules exist and the
## workflow is exactly the one without spike-in: processing.smk / analysis.smk define the plain variants.
##
## Flow:
##   spike_reference      host FASTA/GTF + spike-in FASTA/GTF -> one combined reference (spike-in names prefixed)
##   star_index           (processing.smk) builds on the combined FASTA/GTF when build_star_index is true
##   htseq_count          (processing.smk) counts on the combined GTF, so spike-in genes get counts too
##   spike_idxstats       samtools idxstats of each final BAM
##   spike_summary        host / spike-in read pairs per sample -> spike/spike_counts.tsv (+ per contig)
##   scale_factors        bigWig factors: from spike-in pairs ('tracks' in normalize) or from host pairs only
##   bigwig_sample        the same for per-sample tracks
##   count_matrix         host genes only in the usual matrices; spike-in genes in counts/counts_spike.txt
##   spike_size_factors   DESeq2 size factors (method total) or control genes (method genes)
##   deseq2               the usual DESeq2 run with those size factors ('deseq2' in normalize)

if RUN_PROCESSING:

    ## The combined FASTA is only needed to build the STAR index, so it is temp(): removed once the index
    ## exists (or right away when build_star_index is false and nothing reads it). The GTF is kept: HTSeq reads it.
    rule spike_reference:
        input:
            host_fasta = GENOME_FASTA,
            host_gtf = GTF,
            spike_fasta = SPIKE_SRC_FASTA,
            spike_gtf = SPIKE_SRC_GTF if SPIKE_SRC_GTF else [],
        output:
            fasta = temp(SPIKE_FASTA),
            gtf = SPIKE_GTF,
            contigs = SPIKE_CONTIGS,
            genes = SPIKE_GENES if SPIKE_SRC_GTF else [],
        log:
            OUT + "/logs/spike_reference.log",
        benchmark:
            OUT + "/benchmarks/spike_reference.tsv",
        threads: 1
        params:
            script = SCRIPTS + "/build_spike_reference.py",
            prefix = SPIKE_PREFIX,
            spike_gtf = f"--spike-gtf {SPIKE_SRC_GTF}" if SPIKE_SRC_GTF else "",
            out_genes = f"--out-genes {SPIKE_GENES}" if SPIKE_SRC_GTF else "",
        shell:
            """
            {PY} {params.script} --host-fasta {input.host_fasta} --host-gtf {input.host_gtf} \
                --spike-fasta {input.spike_fasta} {params.spike_gtf} --prefix {params.prefix} \
                --out-fasta {output.fasta} --out-gtf {output.gtf} --out-contigs {output.contigs} \
                {params.out_genes} > {log} 2>&1
            """

    ## Mapped records per contig of the final (MAPQ-filtered, deduplicated) BAM.
    rule spike_idxstats:
        input:
            bam = BAM,
            bai = BAI,
        output:
            stats = SAMPLE_IDXSTATS,
        log:
            OUT + "/{sample}/logs/spike_idxstats.log",
        benchmark:
            OUT + "/{sample}/benchmarks/spike_idxstats.tsv",
        threads: 1
        conda: ENVS + "/ttpipe_tracks.yml"
        shell:
            """
            samtools idxstats {input.bam} > {output.stats} 2> {log}
            """

    rule spike_summary:
        input:
            stats = expand(SAMPLE_IDXSTATS, sample=SAMPLES),
        output:
            counts = SPIKE_COUNTS,
            contigs = SPIKE_CONTIG_COUNTS,
        log:
            OUT + "/logs/spike_summary.log",
        benchmark:
            OUT + "/benchmarks/spike_summary.tsv",
        threads: 1
        params:
            script = SCRIPTS + "/spike_summary.py",
            samples = " ".join(f"--sample {s} {CONDITION_OF[s]} {SAMPLE_IDXSTATS.format(sample=s)}" for s in SAMPLES),
            prefix = SPIKE_PREFIX,
            use = (" ".join(["--use-contigs"] + SPIKE_USE_CONTIGS)) if SPIKE_USE_CONTIGS else "",
        shell:
            """
            {PY} {params.script} {params.samples} --prefix {params.prefix} {params.use} \
                --out-counts {output.counts} --out-contigs {output.contigs} 2> {log}
            """

    ## Condition bigWig factors. 'tracks' in spike_in.normalize: max(spike-in pairs) / this condition's (source
    ## 'spike-in'). Otherwise the usual library-size factor, but on host pairs only (the spike-in reads would
    ## otherwise count as library). config scale_factors still overrides a condition.
    rule scale_factors:
        input:
            counts = SPIKE_COUNTS,
        output:
            table = SCALE_FACTORS,
        log:
            OUT + "/logs/scale_factors.log",
        benchmark:
            OUT + "/benchmarks/scale_factors.tsv",
        threads: 1
        params:
            script = SCRIPTS + "/scale_factors.py",
            basis = "spike" if SPIKE_TRACKS else "host",
            order = " ".join(CONDITIONS),
            overrides = " ".join(f"{k}={v}" for k, v in (config.get("scale_factors") or {}).items()),
        shell:
            """
            {PY} {params.script} --spike-counts {input.counts} --group-by condition --basis {params.basis} \
                --order {params.order} --override {params.overrides} --out {output.table} 2> {log}
            """

    ## Per-sample tracks: the same basis as the condition tracks, one factor per sample (overrides not used).
    rule bigwig_sample:
        input:
            bam = lambda w: BAM.format(sample=w.sample) if w.suffix == "" else STRAND_BAM.format(sample=w.sample, strand=w.suffix[1:]),
            bai = lambda w: BAI.format(sample=w.sample) if w.suffix == "" else STRAND_BAI.format(sample=w.sample, strand=w.suffix[1:]),
            counts = SPIKE_COUNTS,
        output:
            bw = OUT + "/tracks/samples/{sample}{suffix}.bw",
        wildcard_constraints:
            suffix = "|_forward|_reverse",
        log:
            OUT + "/{sample}/logs/bigwig_sample{suffix}.log",
        benchmark:
            OUT + "/{sample}/benchmarks/bigwig_sample{suffix}.tsv",
        threads: 8
        conda: ENVS + "/ttpipe_tracks.yml"
        params:
            bin = config.get("bigwig_bin_size", 50),
            extra = config.get("bigwig_extra", ""),
            script = SCRIPTS + "/scale_factors.py",
            basis = "spike" if SPIKE_TRACKS else "host",
        shell:
            """
            SF=$({PY} {params.script} --spike-counts {input.counts} --group-by sample --basis {params.basis} \
                --factor-of {wildcards.sample} 2> {log})
            bamCoverage --scaleFactor $SF -p {threads} -b {input.bam} -o {output.bw} \
                --binSize {params.bin} {params.extra} >> {log} 2>&1
            """

    ## Host genes only in counts_forward / counts_reverse / counts_deseq2 (same format as without spike-in);
    ## spike-in genes (forward + reverse) in counts_spike.txt when a spike-in GTF is given.
    rule count_matrix:
        input:
            fwd = expand(HTSEQ, sample=SAMPLES, strand="forward"),
            rev = expand(HTSEQ, sample=SAMPLES, strand="reverse"),
            gene_info = GENE_INFO,
        output:
            fwd = COUNTS_FORWARD,
            rev = COUNTS_REVERSE,
            deseq2 = COUNTS_DESEQ2,
            spike = COUNTS_SPIKE if SPIKE_SRC_GTF else [],
        log:
            OUT + "/logs/count_matrix.log",
        benchmark:
            OUT + "/benchmarks/count_matrix.tsv",
        threads: 1
        conda: ENVS + "/ttpipe_count.yml"
        params:
            names = " ".join(SAMPLES),
            mode = COUNT_STRAND,
            script = SCRIPTS + "/make_count_matrix.py",
            prefix = SPIKE_PREFIX,
            spike_out = f"--spike-out {COUNTS_SPIKE}" if SPIKE_SRC_GTF else "",
        shell:
            """
            {PY} {params.script} --names {params.names} --forward {input.fwd} --reverse {input.rev} \
                --count-strand {params.mode} --gene-info {input.gene_info} \
                --out-forward {output.fwd} --out-reverse {output.rev} --out-deseq2 {output.deseq2} \
                --spike-prefix {params.prefix} {params.spike_out} 2> {log}
            """


if SPIKE_DESEQ2:

    rule spike_size_factors:
        input:
            counts = SPIKE_COUNTS if SPIKE_METHOD == "total" else [],
            genes = SPIKE_GENES if SPIKE_METHOD == "genes" else [],
            counts_spike = COUNTS_SPIKE if SPIKE_METHOD == "genes" else [],
        output:
            out = SPIKE_SIZE_OUT,
        log:
            OUT + "/logs/spike_size_factors.log",
        benchmark:
            OUT + "/benchmarks/spike_size_factors.tsv",
        threads: 1
        params:
            script = SCRIPTS + "/spike_size_factors.py",
            args = (f"--spike-counts {SPIKE_COUNTS}" if SPIKE_METHOD == "total" else
                    f"--spike-genes {SPIKE_GENES} --counts-spike {COUNTS_SPIKE} --prefix {SPIKE_PREFIX}"
                    + ((" --use-contigs " + " ".join(SPIKE_USE_CONTIGS)) if SPIKE_USE_CONTIGS else "")),
            method = SPIKE_METHOD,
        shell:
            """
            {PY} {params.script} --method {params.method} {params.args} --out {output.out} 2> {log}
            """

    ## The DESeq2 rule of analysis.smk with the spike-in size factors (it is defined there only when SPIKE_DESEQ2
    ## is false). Method total: sizeFactors() from size_factors.tsv. Method genes: estimateSizeFactors() with the
    ## spike-in genes as controlGenes. The factors used are written to deseq2_size_factors.tsv.
    rule deseq2:
        input:
            counts = COUNTS_MATRIX,
            samples = SAMPLES_DESEQ2,
            sf = SPIKE_SIZE_OUT,
            counts_spike = COUNTS_SPIKE if SPIKE_METHOD == "genes" else [],
        output:
            degs = expand(DEGS, comparison=COMPARISON_NAMES),
            summary = DEG_SUMMARY,
            norm = NORM_COUNTS,
            vst = VST_COUNTS,
            pca = PCA_PDF,
            sf_out = DESEQ2_SIZE_FACTORS,
        params:
            degs_pattern = lambda wildcards: DEGS,
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
            spike = (f"--size-factors {SPIKE_SIZE_FACTORS}" if SPIKE_METHOD == "total" else
                     f"--spike-counts {COUNTS_SPIKE} --spike-method genes --control-genes {SPIKE_CONTROL_GENES}"),
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
                --norm {output.norm} --vst {output.vst} --pca {output.pca} \
                {params.spike} --size-factors-out {output.sf_out} > {log} 2>&1
            """
