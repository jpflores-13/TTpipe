## Shared helpers for TTpipe.snakefile and the rule files.
## Rules carry no cpu/memory/time of their own beyond `threads`; resources live in
## profiles/<site>/config.yaml (set-threads / set-resources), as in MicroCPipe / dietJuicer.
##
## This file is the single source of truth for sample metadata, output paths and targets.
## Every rule file builds its input/output strings from the templates below.

import os
import re
import sys
import pandas as pd

##### Switches #####
RUN_PROCESSING = bool(config.get("run_processing", True))
RUN_ANALYSIS = bool(config.get("run_analysis", True))
RUN_REPORT = bool(config.get("qc_report", True))
if not (RUN_PROCESSING or RUN_ANALYSIS or RUN_REPORT):
    raise ValueError("run_processing, run_analysis and qc_report are all false: nothing to do")

STRANDEDNESS = str(config.get("strandedness", "dUTP"))
if STRANDEDNESS not in ("dUTP", "forward"):
    raise ValueError(f"config 'strandedness' must be 'dUTP' or 'forward', got {STRANDEDNESS!r}")
COUNT_STRAND = str(config.get("count_strand", "sum"))
if COUNT_STRAND not in ("sum", "sense"):
    raise ValueError(f"config 'count_strand' must be 'sum' or 'sense', got {COUNT_STRAND!r}")


def outdir():
    """outputDir, with a relative path placed under $TTPIPE_OUTPUT_ROOT when set."""
    d = os.path.expanduser(str(config.get("outputDir", "output")))
    if not os.path.isabs(d):
        d = os.path.join(os.environ.get("TTPIPE_OUTPUT_ROOT", os.getcwd()), d)
    return d.rstrip("/")


OUT = outdir()
ENVS = os.path.join(workflow.basedir, "..", "envs")
SCRIPTS = os.path.join(workflow.basedir, "..", "scripts")

##### Samplesheet #####
samples = pd.read_csv(config["samplesheet"], sep="\t", dtype=str).fillna("")
_required = {"Name", "Condition", "Replicate"}
if RUN_PROCESSING:
    _required |= {"Read1", "Read2", "Sequencing_Directory"}
_missing = _required - set(samples.columns)
if _missing:
    raise ValueError(f"{config['samplesheet']} is missing column(s): {sorted(_missing)}")
if samples["Name"].duplicated().any():
    raise ValueError(f"{config['samplesheet']}: duplicate Name(s): "
                     f"{sorted(samples['Name'][samples['Name'].duplicated()])}")
for _col in ("Name", "Condition", "Replicate"):
    _bad = [v for v in samples[_col] if not re.fullmatch(r"[A-Za-z0-9._\-]+", v)]
    if _bad:
        raise ValueError(f"{config['samplesheet']}: {_col} values may only contain letters, digits, "
                         f"'.', '_' and '-' (no spaces); offending: {sorted(set(_bad))}")
if RUN_PROCESSING:
    samples["R1"] = [os.path.join(d, f) for d, f in zip(samples["Sequencing_Directory"], samples["Read1"])]
    samples["R2"] = [os.path.join(d, f) for d, f in zip(samples["Sequencing_Directory"], samples["Read2"])]

SAMPLES = list(samples["Name"])
SAMPLE_TABLE = samples.set_index("Name")
CONDITION_OF = SAMPLE_TABLE["Condition"].to_dict()

CONDITIONS = list(dict.fromkeys(samples["Condition"]))
_order = list(config.get("condition_order", []) or [])
if _order:
    if set(_order) != set(CONDITIONS):
        raise ValueError(f"condition_order {_order} must list exactly the conditions in the samplesheet {CONDITIONS}")
    CONDITIONS = _order
SAMPLES_BY_CONDITION = {c: list(samples.loc[samples["Condition"] == c, "Name"]) for c in CONDITIONS}

REFERENCE = str(config.get("reference_condition", CONDITIONS[0]))
if REFERENCE not in CONDITIONS:
    raise ValueError(f"reference_condition {REFERENCE!r} is not a Condition in the samplesheet {CONDITIONS}")
_cmp = config.get("comparisons") or [[c, REFERENCE] for c in CONDITIONS if c != REFERENCE]
COMPARISONS = []  # list of (treatment, reference)
for _pair in _cmp:
    if len(_pair) != 2 or any(p not in CONDITIONS for p in _pair):
        raise ValueError(f"comparisons entries must be [treatment, reference] using samplesheet Conditions; got {_pair}")
    COMPARISONS.append((str(_pair[0]), str(_pair[1])))
COMPARISON_NAMES = [f"{t}_vs_{r}" for t, r in COMPARISONS]
GSEA_COLLECTIONS = dict(config.get("gsea_collections", {}) or {})
STRANDS = ["forward", "reverse"]

wildcard_constraints:
    sample = "|".join(map(re.escape, SAMPLES)),
    condition = "|".join(map(re.escape, CONDITIONS)),
    comparison = "|".join(map(re.escape, COMPARISON_NAMES)) or "NONE",
    collection = "|".join(map(re.escape, GSEA_COLLECTIONS)) or "NONE",
    strand = "forward|reverse",
    mate = "1|2",

##### Reference #####
GTF = config.get("gtf", "")
GENOME_FASTA = config.get("genome_fasta", "")
STAR_INDEX = str(config.get("star_index", "")).rstrip("/")
BUILD_STAR_INDEX = bool(config.get("build_star_index", False))
if RUN_PROCESSING and BUILD_STAR_INDEX and not STAR_INDEX:
    STAR_INDEX = OUT + "/reference/star_index"   # nothing configured: build inside the output dir

def _first_seqname(path, fasta=False):
    """First chromosome name in a FASTA (header) or GTF (first non-comment line); None if unreadable."""
    import gzip
    try:
        op = gzip.open if str(path).endswith(".gz") else open
        with op(path, "rt") as fh:
            for line in fh:
                if fasta and line.startswith(">"):
                    return line[1:].split()[0]
                if not fasta and line.strip() and not line.startswith("#"):
                    return line.split("\t")[0]
    except OSError:
        return None
    return None


def check_reference_names():
    """Fail early when GTF and genome use different chromosome naming styles ('chr1' vs '1'), the trap in the
    source guide (Ensembl FASTA + GENCODE GTF). Only the naming STYLE of the first sequence is compared, so it
    is cheap; it cannot catch every mismatch (e.g. different contig sets)."""
    gtf_name = _first_seqname(GTF) if GTF else None
    genome_name = None
    for fname in ("chrName.txt", "chrNameLength.txt"):          # a STAR index has both
        f = os.path.join(STAR_INDEX, fname) if STAR_INDEX else ""
        if f and os.path.exists(f):
            with open(f) as fh:
                first = fh.readline().split("\t")[0].strip()
            genome_name = first or None
            break
    if genome_name is None and GENOME_FASTA and os.path.exists(GENOME_FASTA):
        genome_name = _first_seqname(GENOME_FASTA, fasta=True)
    if gtf_name and genome_name and gtf_name.startswith("chr") != genome_name.startswith("chr"):
        raise ValueError(f"chromosome names differ between the GTF (first sequence {gtf_name!r}) and the genome "
                         f"(first sequence {genome_name!r}). Use a FASTA/STAR index and a GTF from the same source "
                         "(e.g. both GENCODE, or both UCSC-style 'chr' names). See docs: Setup.")


if RUN_PROCESSING:
    check_reference_names()

##### Spike-in normalization (OPTIONAL, OFF by default, EXPERIMENTAL) #####
## With `spike_in.enabled: false` (the default) nothing below changes the workflow: the same rules, jobs, commands and
## outputs as without this block (tests/test_spike_in.py compares the dry-run job tables). Everything here is
## validated at parse time so that a mistake stops the run before any job starts.
_SPIKE = config.get("spike_in") or {}
SPIKE_ENABLED = bool(_SPIKE.get("enabled", False))
SPIKE_SRC_FASTA = str(_SPIKE.get("fasta") or "")
SPIKE_SRC_GTF = str(_SPIKE.get("gtf") or "")
SPIKE_PREFIX = str(_SPIKE.get("prefix", "spike_") if _SPIKE.get("prefix", "spike_") is not None else "")
SPIKE_METHOD = str(_SPIKE.get("size_factor_method", "total") or "total")
SPIKE_MIN_PAIRS = _SPIKE.get("min_pairs", 10000)


def _as_str_list(v):
    if v is None or v == "":
        return []
    return [str(x) for x in v] if isinstance(v, (list, tuple)) else [str(v)]


SPIKE_USE_CONTIGS = _as_str_list(_SPIKE.get("use_contigs", []))
SPIKE_NORMALIZE = _as_str_list(_SPIKE.get("normalize", ["tracks", "deseq2"]))


def _fasta_names(path):
    """All sequence names (first word of each header) of a FASTA, plain or .gz; None if it cannot be read."""
    import gzip
    try:
        op = gzip.open if str(path).endswith(".gz") else open
        with op(path, "rt") as fh:
            return [ln[1:].split()[0] for ln in fh if ln.startswith(">") and ln[1:].split()]
    except OSError:
        return None


def check_spike_in():
    """Parse-time checks of the spike_in block (only called when it is enabled)."""
    def bad(msg):
        raise ValueError("spike_in: " + msg)
    if not RUN_PROCESSING:
        bad("enabled needs run_processing: true (the spike-in read counts come from the BAMs of this run; "
            "an analysis_counts matrix has none)")
    if not SPIKE_SRC_FASTA:
        bad("enabled is true but 'fasta' is empty: set spike_in.fasta to the spike-in sequences (FASTA, plain or .gz)")
    if not os.path.exists(SPIKE_SRC_FASTA):
        bad(f"fasta {SPIKE_SRC_FASTA!r} does not exist")
    if SPIKE_SRC_GTF and not os.path.exists(SPIKE_SRC_GTF):
        bad(f"gtf {SPIKE_SRC_GTF!r} does not exist")
    if not SPIKE_PREFIX or not re.fullmatch(r"[A-Za-z0-9._\-]+", SPIKE_PREFIX):
        bad(f"prefix {SPIKE_PREFIX!r} must be non-empty and contain only letters, digits, '.', '_' and '-'")
    if SPIKE_METHOD not in ("total", "genes"):
        bad(f"size_factor_method must be 'total' or 'genes', got {SPIKE_METHOD!r}")
    unknown = [x for x in SPIKE_NORMALIZE if x not in ("tracks", "deseq2")]
    if unknown:
        bad(f"normalize may contain 'tracks' and 'deseq2' (or be empty to only count and report); got unknown {unknown}")
    if SPIKE_METHOD == "genes" and "deseq2" in SPIKE_NORMALIZE and not SPIKE_SRC_GTF:
        bad("size_factor_method 'genes' needs spike_in.gtf (the spike-in genes are the control genes of DESeq2)")
    try:
        float(SPIKE_MIN_PAIRS)
        ok = float(SPIKE_MIN_PAIRS) >= 0
    except (TypeError, ValueError):
        ok = False
    if not ok:
        bad(f"min_pairs must be a number >= 0, got {SPIKE_MIN_PAIRS!r}")
    for key, val in (("genome_fasta", GENOME_FASTA), ("gtf", GTF)):
        if not val or not os.path.exists(val):
            bad(f"enabled needs the host {key} ({val!r} is not an existing file): the combined host + spike-in "
                "reference (and the combined GTF used for counting) is built from it")
    names = _fasta_names(SPIKE_SRC_FASTA)
    if names is not None:
        if not names:
            bad(f"no sequences found in fasta {SPIKE_SRC_FASTA!r}")
        missing = [c for c in SPIKE_USE_CONTIGS if c not in names]
        if missing:
            bad(f"use_contigs {missing[:10]} are not sequence names in {SPIKE_SRC_FASTA} (names as in the FASTA, "
                f"without the prefix; it has {len(names)}, e.g. {names[:3]})")
    elif SPIKE_USE_CONTIGS:
        print(f"spike_in: could not read {SPIKE_SRC_FASTA} to check use_contigs", file=sys.stderr)
    if not BUILD_STAR_INDEX:
        idx_names = None
        for fname in ("chrName.txt", "chrNameLength.txt"):
            f = os.path.join(STAR_INDEX, fname) if STAR_INDEX else ""
            if f and os.path.exists(f):
                with open(f) as fh:
                    idx_names = {ln.split("\t")[0].strip() for ln in fh if ln.strip()}
                break
        if idx_names is not None and names:
            lacking = [SPIKE_PREFIX + c for c in names if SPIKE_PREFIX + c not in idx_names]
            if lacking:
                bad(f"the STAR index {STAR_INDEX!r} does not contain the spike-in contigs (e.g. {lacking[:3]}) and "
                    "build_star_index is false. Set build_star_index: true and a NEW star_index directory so that the "
                    "index is built once from the combined host + spike-in reference, or point star_index at an index "
                    "that already contains the prefixed spike-in contigs.")
    elif STAR_INDEX and names:
        for fname in ("chrName.txt", "chrNameLength.txt"):
            f = os.path.join(STAR_INDEX, fname)
            if os.path.exists(f):
                with open(f) as fh:
                    have = {ln.split("\t")[0].strip() for ln in fh if ln.strip()}
                if not {SPIKE_PREFIX + c for c in names} <= have:
                    print(f"spike_in: NOTE: {STAR_INDEX} already holds a STAR index without the spike-in contigs; "
                          "it will be rebuilt (overwritten) from the combined reference. Use a new star_index "
                          "directory to keep the old index.", file=sys.stderr)
                break


if SPIKE_ENABLED:
    check_spike_in()
SPIKE_TRACKS = SPIKE_ENABLED and "tracks" in SPIKE_NORMALIZE            # bigWig scale factors from the spike-in
SPIKE_DESEQ2 = SPIKE_ENABLED and "deseq2" in SPIKE_NORMALIZE and RUN_ANALYSIS   # DESeq2 size factors from the spike-in

##### Output path templates (use with expand() or .format()) #####
## Spike-in (only used when spike_in.enabled is true)
SPIKE_FASTA  = OUT + "/reference/host_spike.fa"           # host FASTA + spike-in FASTA (headers prefixed)
SPIKE_GTF    = OUT + "/reference/host_spike.gtf"          # host GTF + spike-in GTF (seqname, gene_id, transcript_id prefixed; gene lines synthesised if missing)
SPIKE_CONTIGS = OUT + "/reference/spike_contigs.tsv"      # prefixed spike-in contig, length
SPIKE_GENES  = OUT + "/reference/spike_genes.tsv"         # only with spike_in.gtf: gene_id, contig, start, end, strand
SAMPLE_IDXSTATS = OUT + "/{sample}/stats/{sample}_idxstats.txt"      # samtools idxstats of the final BAM
SPIKE_COUNTS = OUT + "/spike/spike_counts.tsv"            # sample, condition, host_pairs, spike_pairs_all, spike_pairs_used, spike_pct
SPIKE_CONTIG_COUNTS = OUT + "/spike/spike_contig_counts.tsv"   # spike-in contig x sample read pairs
SPIKE_SIZE_FACTORS = OUT + "/spike/size_factors.tsv"      # method total: sample, size_factor (for DESeq2)
SPIKE_CONTROL_GENES = OUT + "/spike/control_genes.txt"    # method genes: control gene ids for DESeq2
SPIKE_SIZE_OUT = SPIKE_SIZE_FACTORS if SPIKE_METHOD == "total" else SPIKE_CONTROL_GENES
COUNTS_SPIKE = OUT + "/counts/counts_spike.txt"           # spike-in genes x samples (only with spike_in.gtf)
## What the reference-dependent rules read: the plain host files, or the combined host + spike-in ones
STAR_FASTA = SPIKE_FASTA if SPIKE_ENABLED else GENOME_FASTA
STAR_GTF = SPIKE_GTF if SPIKE_ENABLED else GTF
HTSEQ_GTF = SPIKE_GTF if SPIKE_ENABLED else GTF
## Per sample
FASTQC_ZIP   = OUT + "/qc/fastqc/{sample}_R{mate}_fastqc.zip"
TRIM_R       = OUT + "/{sample}/trim/{sample}_R{mate}.fastq.gz"          # temp, only if trim: true
FASTP_JSON   = OUT + "/{sample}/stats/{sample}_fastp.json"
STAR_BAM     = OUT + "/{sample}/star/{sample}_Aligned.sortedByCoord.out.bam"   # temp
STAR_LOG     = OUT + "/{sample}/stats/{sample}_Log.final.out"
STAR_GENECNT = OUT + "/{sample}/stats/{sample}_ReadsPerGene.out.tab"
MAPQ_BAM     = OUT + "/{sample}/{sample}_mapq.bam"                        # temp
MAPQ_STATS   = OUT + "/{sample}/stats/{sample}_mapq.txt"                  # one line: "<sample>\t<read pairs after MAPQ filter>"
DEDUP_METRICS = OUT + "/{sample}/stats/{sample}_dedupMetrics.txt"
BAM          = OUT + "/{sample}/{sample}.bam"                             # final: MAPQ-filtered, deduplicated, indexed
BAI          = BAM + ".bai"
STRAND_BAM   = OUT + "/{sample}/{sample}_{strand}.bam"                    # forward = RNA transcribed from the + genomic strand, reverse = from the - strand (NOT "sense": a + gene's sense reads are in forward, a - gene's in reverse)
STRAND_BAI   = STRAND_BAM + ".bai"
FINAL_STATS  = OUT + "/{sample}/stats/{sample}_final_pairs.txt"           # one line: "<sample>\t<read pairs in final BAM>"
HTSEQ        = OUT + "/{sample}/counts/{sample}_{strand}.counts.txt"      # htseq-count output, gene_id<TAB>count, with __ summary rows
## Per condition (replicates merged)
MERGED_BAM   = OUT + "/merged/{condition}.bam"
MERGED_STRAND_BAM = OUT + "/merged/{condition}_{strand}.bam"
SCALE_FACTORS = OUT + "/tracks/scale_factors.tsv"                         # columns: condition, mapped_pairs, scale_factor, source
BIGWIG       = OUT + "/tracks/{condition}.bw"
BIGWIG_STRAND = OUT + "/tracks/{condition}_{strand}.bw"
SAMPLE_BIGWIG = OUT + "/tracks/samples/{sample}.bw"                       # only if sample_tracks: true
SAMPLE_BIGWIG_STRAND = OUT + "/tracks/samples/{sample}_{strand}.bw"
## Counts
GENE_INFO    = OUT + "/counts/gene_info.tsv"          # gene_id, gene_name, gene_type, chr, start, end, strand, length
COUNTS_FORWARD = OUT + "/counts/counts_forward.txt"   # genes x samples, header = sample Names, first col = gene_id (HTSeq __ rows removed)
COUNTS_REVERSE = OUT + "/counts/counts_reverse.txt"
COUNTS_DESEQ2  = OUT + "/counts/counts_deseq2.txt"    # per count_strand: sum (fwd+rev) or sense
SAMPLES_DESEQ2 = OUT + "/counts/samples_deseq2.txt"   # columns: sample, condition, replicate
COUNTS_MATRIX  = config.get("analysis_counts") or COUNTS_DESEQ2   # what the analysis reads
if RUN_ANALYSIS and not RUN_PROCESSING and not config.get("analysis_counts") and not os.path.exists(COUNTS_DESEQ2):
    raise ValueError("run_processing is false and analysis_counts is empty, and "
                     f"{COUNTS_DESEQ2} does not exist: set analysis_counts to a count matrix")
## Analysis
DEG_DIR      = OUT + "/analysis/deseq2"
DEGS         = DEG_DIR + "/degs_{comparison}.csv"     # full DESeq2 results table, rows = genes (gene_id), as.data.frame(results)
DEG_SUMMARY  = DEG_DIR + "/deg_summary.tsv"           # comparison, n_up, n_down, n_tested
NORM_COUNTS  = DEG_DIR + "/normalized_counts.tsv"
VST_COUNTS   = DEG_DIR + "/vst.tsv"
DESEQ2_SIZE_FACTORS = DEG_DIR + "/deseq2_size_factors.tsv"   # only with spike-in size factors: sample, size_factor, method
PCA_PDF      = OUT + "/analysis/figures/pca.pdf"
GSEA_DIR     = OUT + "/analysis/gsea"
RNK          = GSEA_DIR + "/{comparison}.rnk"         # gene_name <TAB> Wald stat, sorted descending
GSEA_RESULT  = GSEA_DIR + "/{collection}_{comparison}.tsv"
SANKEY_PDF   = OUT + "/analysis/figures/sankey.pdf"
SANKEY_TABLE = OUT + "/analysis/deseq2/sankey_categories.tsv"   # gene_id + Up/Down/n.s. per comparison
## QC
SAMPLE_QC    = OUT + "/qc/sample_qc.tsv"              # one row per sample, see scripts/qc_summary.py
QC_REPORT    = OUT + "/qc/TTpipe_qc_report.html"
MULTIQC_REPORT = OUT + "/qc/TTpipe_multiqc.html"
## Results explorer (static interactive HTML: volcano/MA, gene search, tables, GSEA); only if `results_explorer: true`
RESULTS_EXPLORER = OUT + "/analysis/TTpipe_results_explorer.html"
## Stranded signal page (locus tracks, stranded metagene profiles, antisense-fraction table); only if `stranded_plots: true`
STRANDED_SIGNAL = OUT + "/analysis/TTpipe_stranded_signal.html"


##### Targets #####
def processing_targets():
    t = []
    t += expand(FASTQC_ZIP, sample=SAMPLES, mate=["1", "2"])
    t += expand(BAI, sample=SAMPLES)
    t += expand(STRAND_BAI, sample=SAMPLES, strand=STRANDS)
    t += expand(BIGWIG, condition=CONDITIONS)
    t += expand(BIGWIG_STRAND, condition=CONDITIONS, strand=STRANDS)
    if config.get("sample_tracks", False):
        t += expand(SAMPLE_BIGWIG, sample=SAMPLES)
        t += expand(SAMPLE_BIGWIG_STRAND, sample=SAMPLES, strand=STRANDS)
    t += [COUNTS_FORWARD, COUNTS_REVERSE, COUNTS_DESEQ2, GENE_INFO, SCALE_FACTORS]
    if SPIKE_ENABLED:
        t += [SPIKE_COUNTS, SPIKE_CONTIG_COUNTS, SPIKE_CONTIGS]
        if SPIKE_SRC_GTF:
            t.append(COUNTS_SPIKE)
    return t


def analysis_targets():
    t = [SAMPLES_DESEQ2, GENE_INFO, DEG_SUMMARY, NORM_COUNTS, VST_COUNTS, PCA_PDF, SANKEY_PDF, SANKEY_TABLE]
    t += expand(DEGS, comparison=COMPARISON_NAMES)
    if SPIKE_DESEQ2:
        t.append(DESEQ2_SIZE_FACTORS)
    if GSEA_COLLECTIONS:
        t += expand(GSEA_RESULT, collection=list(GSEA_COLLECTIONS), comparison=COMPARISON_NAMES)
    if config.get("results_explorer", True):
        t.append(RESULTS_EXPLORER)
    return t


def stranded_targets():
    """The stranded signal page, when `stranded_plots` is true (default) AND the stranded bigWigs are made in this run
    or already exist (an analysis-only run on an output folder that was processed earlier). Without bigWigs the page is
    skipped quietly: an analysis-only run on a count matrix has nothing to draw."""
    if not config.get("stranded_plots", True):
        return []
    if RUN_PROCESSING:
        return [STRANDED_SIGNAL]
    if RUN_ANALYSIS:
        needed = (expand(BIGWIG_STRAND, condition=CONDITIONS, strand=STRANDS)
                  + [SCALE_FACTORS, GENE_INFO, COUNTS_FORWARD, COUNTS_REVERSE])
        if all(os.path.exists(f) for f in needed):
            return [STRANDED_SIGNAL]
    return []


def report_targets():
    return [QC_REPORT]


def all_targets():
    t = []
    if RUN_PROCESSING:
        t += processing_targets()
    if RUN_ANALYSIS:
        t += analysis_targets()
    t += stranded_targets()
    if RUN_REPORT:
        t += report_targets()
    if RUN_PROCESSING and config.get("multiqc", False):
        t.append(MULTIQC_REPORT)
    return t


##### Email notifications (see `notify:` in the config) #####
NOTIFY = config.get("notify", {}) or {}
EMAIL = os.environ.get("TTPIPE_EMAIL") or NOTIFY.get("email") or ""


def notify(event, body, run_name):
    """Send one message for a workflow event, if an address and the event are configured."""
    if not EMAIL or event not in NOTIFY.get("events", ["start", "success", "error"]):
        return
    import shlex, socket, subprocess
    subject = f"{run_name}: {event}"
    body = (f"Run:     {run_name}\nOutput:  {OUT}\nHost:    {socket.gethostname()}\n\n{body}\n")
    cmd = NOTIFY.get("command", "mail -s {subject} {email}").format(
        subject=shlex.quote(subject), email=shlex.quote(EMAIL))
    try:
        subprocess.run(cmd, shell=True, input=body, text=True, timeout=60, check=True)
    except Exception as e:  # a failed email must never fail the workflow
        print(f"Could not send '{event}' email to {EMAIL}: {e}", file=sys.stderr)


def log_tail(path, n=40):
    """Last n lines of the Snakemake log, for the failure email."""
    try:
        with open(path) as fh:
            return "".join(fh.readlines()[-n:])
    except (OSError, TypeError):
        return "(Snakemake log unavailable)"
