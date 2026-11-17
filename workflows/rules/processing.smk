## TTpipe processing rules: FASTQ -> QC -> [fastp] -> STAR -> MAPQ filter -> Picard dedup -> strand BAMs
## -> replicate-merged BAMs -> scale factors -> bigWigs; HTSeq counts -> count matrices.
## Logic and parameters follow the TT-seq guide, sections 4.2-4.8 (deliberate differences are commented).
## Paths come from rules/common.smk. Rules set only `threads:`; memory / time live in profiles/<site>/config.yaml.
## Rule names (profiles key on them): fastqc, fastp, star_index, star_align, mapq_filter, mark_duplicates,
## split_strands, merge_bams, scale_factors, bigwig, bigwig_sample, gene_info, htseq_count, count_matrix.

TRIM = bool(config.get("trim", False))
PY = "$CONDA_PREFIX/bin/python"

if RUN_PROCESSING and not STAR_INDEX:
    raise ValueError("config 'star_index' is empty: set it to a STAR genome directory, "
                     "or set build_star_index: true to build one")


def raw_reads(w):
    return [SAMPLE_TABLE.loc[w.sample, "R1"], SAMPLE_TABLE.loc[w.sample, "R2"]]


def star_reads(w):
    """Trimmed reads when trim is true, raw FASTQs otherwise."""
    if TRIM:
        return [TRIM_R.format(sample=w.sample, mate=m) for m in ("1", "2")]
    return raw_reads(w)


## Strand assignment of paired reads (samtools -f/-F), per library type. 'forward' here means RNA
## transcribed from the + genomic strand (NOT "sense"; a - gene's sense reads are in 'reverse').
##   dUTP (guide): read 1 antisense, read 2 sense.
##     forward = read2 not reversed (-f 128 -F 16) + read1 reversed (-f 80)
##     reverse = read2 reversed (-f 144)           + read1 not reversed (-f 64 -F 16)
##   forward-stranded (TTpipe addition; read 1 sense): the two sets swap.
##     forward = read1 not reversed (-f 64 -F 16)  + read2 reversed (-f 144)
##     reverse = read1 reversed (-f 80)            + read2 not reversed (-f 128 -F 16)
STRAND_FLAGS = {
    "dUTP": {"forward": ("-f 128 -F 16", "-f 80"), "reverse": ("-f 144", "-f 64 -F 16")},
    "forward": {"forward": ("-f 64 -F 16", "-f 144"), "reverse": ("-f 80", "-f 128 -F 16")},
}


if RUN_PROCESSING:

    ##### QC #####
    ## FastQC always runs on the raw FASTQs (as the guide). It names outputs after the input file, so the
    ## FASTQ is symlinked as {sample}_R{mate}.fastq.gz in a scratch dir to get FASTQC_ZIP's name.
    rule fastqc:
        input:
            fq = lambda w: SAMPLE_TABLE.loc[w.sample, "R" + w.mate],
        output:
            zip = FASTQC_ZIP,
        log:
            OUT + "/{sample}/logs/fastqc_R{mate}.log",
        benchmark:
            OUT + "/{sample}/benchmarks/fastqc_R{mate}.tsv",
        threads: 4
        conda: ENVS + "/ttpipe_align.yml"
        params:
            outdir = OUT + "/qc/fastqc",
            tmp = OUT + "/qc/fastqc/tmp_{sample}_R{mate}",
        shell:
            """
            rm -rf {params.tmp}; mkdir -p {params.tmp}
            ln -s "$(readlink -f {input.fq})" {params.tmp}/{wildcards.sample}_R{wildcards.mate}.fastq.gz
            fastqc -t {threads} -o {params.outdir} {params.tmp}/{wildcards.sample}_R{wildcards.mate}.fastq.gz > {log} 2>&1
            rm -rf {params.tmp}
            """

    ## Optional adapter / quality trimming (not in the guide; config trim: true). fastp defaults
    ## (automatic PE adapter detection, quality filter, length filter) plus fastp_extra.
    if TRIM:
        rule fastp:
            input:
                reads = raw_reads,
            output:
                r1 = temp(TRIM_R.replace("{mate}", "1")),
                r2 = temp(TRIM_R.replace("{mate}", "2")),
                json = FASTP_JSON,
                html = OUT + "/{sample}/stats/{sample}_fastp.html",
            log:
                OUT + "/{sample}/logs/fastp.log",
            benchmark:
                OUT + "/{sample}/benchmarks/fastp.tsv",
            threads: 8
            conda: ENVS + "/ttpipe_align.yml"
            params:
                extra = config.get("fastp_extra", ""),
            shell:
                """
                fastp -i {input.reads[0]} -I {input.reads[1]} -o {output.r1} -O {output.r2} \
                    -w {threads} -j {output.json} -h {output.html} {params.extra} > {log} 2>&1
                """

    ##### Alignment #####
    ## Only defined when build_star_index is true; otherwise star_align just needs the existing index.
    if BUILD_STAR_INDEX:
        rule star_index:
            input:
                fasta = STAR_FASTA,       # = genome_fasta and gtf, or the combined host + spike-in files (spike_in.enabled)
                gtf = STAR_GTF,
            output:
                sa = STAR_INDEX + "/SA",
                chrlen = STAR_INDEX + "/chrNameLength.txt",
            log:
                OUT + "/logs/star_index.log",
            benchmark:
                OUT + "/benchmarks/star_index.tsv",
            threads: 16
            conda: ENVS + "/ttpipe_align.yml"
            params:
                index = STAR_INDEX,
                overhang = config.get("sjdb_overhang", 100),
                prefix = OUT + "/logs/star_index_",
            shell:
                ## TTpipe change: --sjdbOverhang is passed explicitly (guide: STAR default, which is 100).
                """
                mkdir -p {params.index}
                STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.index} \
                    --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} \
                    --sjdbOverhang {params.overhang} --outFileNamePrefix {params.prefix} > {log} 2>&1
                """

    ## Guide parameters: --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --quantMode GeneCounts.
    ## Housekeeping as the guide: SJ.out.tab, Log.out, Log.progress.out deleted; Log.final.out and
    ## ReadsPerGene.out.tab kept under stats/. TTpipe addition: per-sample --outTmpDir (removed afterwards).
    rule star_align:
        input:
            reads = star_reads,
            index = STAR_INDEX + "/SA",
        output:
            bam = temp(STAR_BAM),
            log_final = STAR_LOG,
            genecnt = STAR_GENECNT,
        log:
            OUT + "/{sample}/logs/star_align.log",
        benchmark:
            OUT + "/{sample}/benchmarks/star_align.tsv",
        threads: 16
        conda: ENVS + "/ttpipe_align.yml"
        params:
            index = STAR_INDEX,
            prefix = OUT + "/{sample}/star/{sample}_",
            tmp = OUT + "/{sample}/star/{sample}_STARtmp",
            extra = config.get("star_extra", ""),
        shell:
            """
            rm -rf {params.tmp}
            STAR --runThreadN {threads} --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat \
                --quantMode GeneCounts --genomeDir {params.index} --outFileNamePrefix {params.prefix} \
                --outTmpDir {params.tmp} {params.extra} \
                --readFilesIn {input.reads[0]} {input.reads[1]} > {log} 2>&1
            mv {params.prefix}Log.final.out {output.log_final}
            mv {params.prefix}ReadsPerGene.out.tab {output.genecnt}
            rm -rf {params.tmp} {params.prefix}STARgenome {params.prefix}STARpass1
            rm -f {params.prefix}SJ.out.tab {params.prefix}Log.out {params.prefix}Log.progress.out
            """

    ## samtools view -h -S -b -q <min_mapq> (guide: 10). Read pairs after the filter = records / 2.
    rule mapq_filter:
        input:
            bam = STAR_BAM,
        output:
            bam = temp(MAPQ_BAM),
            stats = MAPQ_STATS,
        log:
            OUT + "/{sample}/logs/mapq_filter.log",
        benchmark:
            OUT + "/{sample}/benchmarks/mapq_filter.tsv",
        threads: 8
        conda: ENVS + "/ttpipe_align.yml"
        params:
            mapq = config.get("min_mapq", 10),
        shell:
            """
            samtools view -@ {threads} -h -S -b -q {params.mapq} -o {output.bam} {input.bam} 2> {log}
            n=$(samtools view -@ {threads} -c {output.bam})
            printf '%s\\t%s\\n' {wildcards.sample} $((n / 2)) > {output.stats}
            """

    ## Picard MarkDuplicates (guide: REMOVE_DUPLICATES=true) then samtools index.
    ## No Java heap is set here (rules carry no memory); Picard uses the JVM default. Add
    ## e.g. picard_extra: '-Xmx12g' only if your picard wrapper accepts it (not verified).
    rule mark_duplicates:
        input:
            bam = MAPQ_BAM,
        output:
            bam = BAM,
            bai = BAI,
            metrics = DEDUP_METRICS,
            stats = FINAL_STATS,
        log:
            OUT + "/{sample}/logs/mark_duplicates.log",
        benchmark:
            OUT + "/{sample}/benchmarks/mark_duplicates.tsv",
        threads: 2
        conda: ENVS + "/ttpipe_align.yml"
        params:
            dedup_flag = "true" if config.get("remove_duplicates", True) else "false",
            extra = config.get("picard_extra", ""),
            tmp = OUT + "/{sample}/tmp_picard",
        shell:
            """
            rm -rf {params.tmp}; mkdir -p {params.tmp}
            picard MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} \
                REMOVE_DUPLICATES={params.dedup_flag} TMP_DIR={params.tmp} {params.extra} > {log} 2>&1
            rm -rf {params.tmp}
            samtools index -@ {threads} {output.bam}
            n=$(samtools view -@ {threads} -c {output.bam})
            printf '%s\\t%s\\n' {wildcards.sample} $((n / 2)) > {output.stats}
            """

    ## Strand split: the two flag sets in STRAND_FLAGS (see the table above), then merge.
    ## Both inputs of the merge are coordinate-sorted subsets of a sorted BAM, so the merge is already
    ## sorted: the guide's extra `samtools sort` is skipped (same result).
    rule split_strands:
        input:
            bam = BAM,
            bai = BAI,
        output:
            bam = STRAND_BAM,
            bai = STRAND_BAI,
        log:
            OUT + "/{sample}/logs/split_strands_{strand}.log",
        benchmark:
            OUT + "/{sample}/benchmarks/split_strands_{strand}.tsv",
        threads: 8
        conda: ENVS + "/ttpipe_align.yml"
        params:
            flags1 = lambda w: STRAND_FLAGS[STRANDEDNESS][w.strand][0],
            flags2 = lambda w: STRAND_FLAGS[STRANDEDNESS][w.strand][1],
            tmp = OUT + "/{sample}/tmp_split_{strand}",
        shell:
            """
            rm -rf {params.tmp}; mkdir -p {params.tmp}
            trap 'rm -rf {params.tmp}' EXIT
            samtools view -@ {threads} -b {params.flags1} -o {params.tmp}/part1.bam {input.bam} 2> {log}
            samtools view -@ {threads} -b {params.flags2} -o {params.tmp}/part2.bam {input.bam} 2>> {log}
            samtools merge -@ {threads} -f {output.bam} {params.tmp}/part1.bam {params.tmp}/part2.bam 2>> {log}
            samtools index -@ {threads} {output.bam} 2>> {log}
            """

    ##### Tracks #####
    ## Replicates of a condition merged: all reads (suffix '') and per strand (suffix '_forward' / '_reverse').
    ## Output names equal MERGED_BAM / MERGED_STRAND_BAM; one rule serves both. Merged BAMs are for tracks only.
    ## Inputs are sorted, so the merge is sorted (guide sorts again; same result).
    rule merge_bams:
        input:
            bams = lambda w: [(BAM.format(sample=s) if w.suffix == "" else
                               STRAND_BAM.format(sample=s, strand=w.suffix[1:]))
                              for s in SAMPLES_BY_CONDITION[w.condition]],
        output:
            bam = OUT + "/merged/{condition}{suffix}.bam",
            bai = OUT + "/merged/{condition}{suffix}.bam.bai",
        wildcard_constraints:
            suffix = "|_forward|_reverse",
        log:
            OUT + "/logs/merge_bams_{condition}{suffix}.log",
        benchmark:
            OUT + "/benchmarks/merge_bams_{condition}{suffix}.tsv",
        threads: 8
        conda: ENVS + "/ttpipe_tracks.yml"
        shell:
            """
            if [ $(echo {input.bams} | wc -w) -eq 1 ]; then
                cp {input.bams} {output.bam} 2> {log}
            else
                samtools merge -@ {threads} -f {output.bam} {input.bams} 2> {log}
            fi
            samtools index -@ {threads} {output.bam} 2>> {log}
            """

    ## Scale factor per condition = (largest mapped read-pair count) / (this one's), on the merged all-reads
    ## BAMs, pairs = records / 2 (guide). config scale_factors overrides a condition (source 'config').
    ## The maximum is taken over all conditions, including overridden ones.
    if not SPIKE_ENABLED:   # the spike-in variants of this rule are in rules/spike.smk
        rule scale_factors:
            input:
                bams = expand(OUT + "/merged/{condition}.bam", condition=CONDITIONS),
            output:
                table = SCALE_FACTORS,
            log:
                OUT + "/logs/scale_factors.log",
            benchmark:
                OUT + "/benchmarks/scale_factors.tsv",
            threads: 1
            conda: ENVS + "/ttpipe_tracks.yml"
            params:
                cond_bams = " ".join(f"{c}={OUT}/merged/{c}.bam" for c in CONDITIONS),
                overrides = " ".join(f"{k}={v}" for k, v in (config.get("scale_factors") or {}).items()),
                script = SCRIPTS + "/scale_factors.py",
                tmp = OUT + "/tracks/scale_factors_pairs.tmp",
            shell:
                """
            : > {params.tmp}
            for kv in {params.cond_bams}; do
                c=${{kv%%=*}}; b=${{kv#*=}}
                n=$(samtools view -@ {threads} -c "$b")
                printf '%s\\t%s\\n' "$c" $((n / 2)) >> {params.tmp}
            done
            {PY} {params.script} --pairs {params.tmp} --override {params.overrides} --out {output.table} 2> {log}
            rm -f {params.tmp}
            """

    ## bamCoverage --scaleFactor SF -p threads -b BAM -o BW (guide; no other normalization).
    ## --binSize is passed explicitly (default 50 = deepTools' own default, so the same as the guide).
    rule bigwig:
        input:
            bam = OUT + "/merged/{condition}{suffix}.bam",
            bai = OUT + "/merged/{condition}{suffix}.bam.bai",
            sf = SCALE_FACTORS,
        output:
            bw = OUT + "/tracks/{condition}{suffix}.bw",
        wildcard_constraints:
            suffix = "|_forward|_reverse",
        log:
            OUT + "/logs/bigwig_{condition}{suffix}.log",
        benchmark:
            OUT + "/benchmarks/bigwig_{condition}{suffix}.tsv",
        threads: 8
        conda: ENVS + "/ttpipe_tracks.yml"
        params:
            bin = config.get("bigwig_bin_size", 50),
            extra = config.get("bigwig_extra", ""),
        shell:
            """
            SF=$(awk -F'\\t' -v c={wildcards.condition} '$1==c {{print $3}}' {input.sf})
            [ -n "$SF" ] || {{ echo "no scale factor for {wildcards.condition}" > {log}; exit 1; }}
            bamCoverage --scaleFactor $SF -p {threads} -b {input.bam} -o {output.bw} \
                --binSize {params.bin} {params.extra} > {log} 2>&1
            """

    ## Per-sample tracks (config sample_tracks: true). TTpipe addition, not in the guide.
    ## A condition's scale factor does not apply to a single replicate, so the factor here is
    ## (largest final read-pair count over ALL samples) / (this sample's), from the *_final_pairs.txt files,
    ## the same max/this rule as the condition tracks. config scale_factors is not used for these.
    if not SPIKE_ENABLED:   # the spike-in variants of this rule are in rules/spike.smk
        rule bigwig_sample:
            input:
                bam = lambda w: BAM.format(sample=w.sample) if w.suffix == "" else STRAND_BAM.format(sample=w.sample, strand=w.suffix[1:]),
                bai = lambda w: BAI.format(sample=w.sample) if w.suffix == "" else STRAND_BAI.format(sample=w.sample, strand=w.suffix[1:]),
                stats = expand(FINAL_STATS, sample=SAMPLES),
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
            shell:
                """
            SF=$({PY} {params.script} --pairs {input.stats} --factor-of {wildcards.sample} 2> {log})
            bamCoverage --scaleFactor $SF -p {threads} -b {input.bam} -o {output.bw} \
                --binSize {params.bin} {params.extra} >> {log} 2>&1
            """

    ##### Counts #####
    ## Guide: htseq-count -f bam --stranded no --nonunique all --order pos --type gene on each strand BAM.
    ## --order pos with paired-end data buffers reads whose mate has not been seen yet in memory (this is
    ## what HTSeq's --max-reads-in-buffer limits); both mates of a pair are in the same strand BAM.
    ## The guide gives all samples to one htseq-count call; here it is one call per sample and strand
    ## (same counts, parallelizable), merged in count_matrix.
    rule htseq_count:
        input:
            bam = STRAND_BAM,
            gtf = HTSEQ_GTF,          # = gtf, or the combined host + spike-in GTF (spike_in.enabled)
        output:
            counts = HTSEQ,
        log:
            OUT + "/{sample}/logs/htseq_count_{strand}.log",
        benchmark:
            OUT + "/{sample}/benchmarks/htseq_count_{strand}.tsv",
        threads: 1
        conda: ENVS + "/ttpipe_count.yml"
        params:
            nonunique = config.get("htseq_nonunique", "all"),
            feature = config.get("htseq_feature", "gene"),
            idattr = config.get("htseq_id_attribute", "gene_id"),
        shell:
            """
            htseq-count -f bam --stranded no --nonunique {params.nonunique} --order pos \
                --type {params.feature} --idattr {params.idattr} {input.bam} {input.gtf} \
                > {output.counts} 2> {log}
            """

    ## Per-sample HTSeq files -> gene x sample matrices (summary rows removed). counts_deseq2 per count_strand.
    if not SPIKE_ENABLED:   # the spike-in variants of this rule are in rules/spike.smk
        rule count_matrix:
            input:
                fwd = expand(HTSEQ, sample=SAMPLES, strand="forward"),
                rev = expand(HTSEQ, sample=SAMPLES, strand="reverse"),
                gene_info = GENE_INFO,
            output:
                fwd = COUNTS_FORWARD,
                rev = COUNTS_REVERSE,
                deseq2 = COUNTS_DESEQ2,
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
            shell:
                """
            {PY} {params.script} --names {params.names} --forward {input.fwd} --reverse {input.rev} \
                --count-strand {params.mode} --gene-info {input.gene_info} \
                --out-forward {output.fwd} --out-reverse {output.rev} --out-deseq2 {output.deseq2} 2> {log}
            """


## Gene table from the GTF `gene` lines. Defined even when run_processing is false (analysis needs it).
rule gene_info:
    input:
        gtf = GTF,
    output:
        table = GENE_INFO,
    log:
        OUT + "/logs/gene_info.log",
    benchmark:
        OUT + "/benchmarks/gene_info.tsv",
    threads: 1
    conda: ENVS + "/ttpipe_count.yml"
    params:
        idattr = config.get("htseq_id_attribute", "gene_id"),
        script = SCRIPTS + "/gene_info.py",
    shell:
        """
        {PY} {params.script} --gtf {input.gtf} --out {output.table} --id-attribute {params.idattr} 2> {log}
        """
