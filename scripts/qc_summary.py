#!/usr/bin/env python3
"""Per-sample QC table for TTpipe (one row per sample) -> <outputDir>/qc/sample_qc.tsv.

Usage:
    python scripts/qc_summary.py --outdir OUTDIR --samplesheet samplesheet.txt \
        --config config/TTpipeConfig.yaml [more.yaml ...] [--counts COUNT_MATRIX] --out OUTDIR/qc/sample_qc.tsv

File discovery. The script rediscovers the files of each sample by the SAME naming scheme as
workflows/rules/common.smk (OUT = --outdir, <s> = sample Name from the samplesheet):
    OUT/<s>/stats/<s>_Log.final.out                STAR_LOG
    OUT/<s>/stats/<s>_ReadsPerGene.out.tab         STAR_GENECNT
    OUT/<s>/stats/<s>_mapq.txt                     MAPQ_STATS    "<sample>\\t<read pairs after MAPQ filter>"
    OUT/<s>/stats/<s>_dedupMetrics.txt             DEDUP_METRICS (Picard MarkDuplicates metrics file)
    OUT/<s>/stats/<s>_final_pairs.txt              FINAL_STATS   "<sample>\\t<read pairs in the final BAM>"
    OUT/counts/counts_forward.txt, counts_reverse.txt   COUNTS_FORWARD / COUNTS_REVERSE (genes x samples)
    OUT/counts/gene_info.tsv                       GENE_INFO (gene_id, gene_type, strand, ...)
    --counts (config `analysis_counts`, else OUT/counts/counts_deseq2.txt)   the matrix used for replicate correlation
Every input is optional. A metric whose input is missing is left BLANK in the output (never
estimated), so processing-only, analysis-only (only a count matrix) and full runs all work.

Columns written (see METRIC_COLUMNS for the order):
    sample, condition, replicate
    input_pairs, input_pairs_million   STAR "Number of input reads" (read pairs for paired-end data)
    unique_pct                         STAR "Uniquely mapped reads %"
    multi_pct                          STAR "% mapped to multiple loci" + "% mapped to too many loci"
    too_short_pct                      STAR "% of reads unmapped: too short"
    other_unmapped_pct                 STAR "% unmapped: too many mismatches" + "% unmapped: other"
    mapq_pairs, final_pairs            read pairs after the MAPQ filter / in the final deduplicated BAM
    final_pct_of_input                 final_pairs / input_pairs * 100
    duplicate_pct                      Picard PERCENT_DUPLICATION * 100 (as reported by Picard)
    sense_pct                          % of counted reads on the annotated gene strand; HTSeq-based when the
                                       strand-split matrices + gene_info exist, else STAR-based (sense_source)
    sense_pct_htseq, sense_pct_star    the two estimates separately
    rrna_pct                           % of counted reads in genes of type rRNA / Mt_rRNA (GENE_INFO gene_type)
    counted_reads                      total reads in the count matrix used for rrna_pct
    replicate_correlation              mean Pearson r of log2(CPM+1) vs the other replicates of the same condition

Definitions and caveats
  * sense_pct (HTSeq): forward matrix = RNA transcribed from the + genomic strand, reverse = from the - strand
    (common.smk). For a + gene the sense reads are in forward, for a - gene in reverse. Genes without a +/-
    strand are ignored. Overlapping genes on opposite strands make this an approximation.
  * sense_pct (STAR): ReadsPerGene columns are unstranded, "read 1 same strand as gene", "read 1 opposite
    strand". Config strandedness 'forward' => column 3 is sense; 'dUTP' => column 4 is sense. STAR counts only
    exonic, unambiguous reads, which are a small part of nascent TT-seq signal, so the HTSeq (gene-body)
    estimate is preferred when available.
  * rrna_pct: only genes annotated rRNA / Mt_rRNA in the GTF count. Mouse (mm10) GENCODE contains very few
    rRNA loci (most rRNA repeats are not annotated), so this UNDERESTIMATES true rRNA contamination.
  * Picard's metrics-file layout (a '## METRICS CLASS' line, a header line, then one row per library) is parsed
    by header name. The exact meaning of PERCENT_DUPLICATION for paired-end data is Picard's; it is used as is.
  * STAR Log.final.out labels are matched on their text (e.g. 'Uniquely mapped reads %'); they were written from
    memory of STAR 2.7 output and tested only on synthetic files, so check them against a real log.
"""
from __future__ import annotations

import argparse
import math
import os
import re
import sys
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import yaml

RRNA_TYPES = {"rRNA", "Mt_rRNA"}

METRIC_COLUMNS = [
    "sample", "condition", "replicate", "input_pairs", "input_pairs_million", "unique_pct", "multi_pct",
    "too_short_pct", "other_unmapped_pct", "mapq_pairs", "final_pairs", "final_pct_of_input", "duplicate_pct",
    "sense_pct", "sense_source", "sense_pct_htseq", "sense_pct_star", "rrna_pct", "counted_reads",
    "replicate_correlation",
]


# --------------------------------------------------------------------------
# Paths (mirror common.smk)
# --------------------------------------------------------------------------

def sample_paths(outdir: str, s: str) -> Dict[str, str]:
    st = os.path.join(outdir, s, "stats")
    return {
        "star_log": os.path.join(st, f"{s}_Log.final.out"),
        "star_genecnt": os.path.join(st, f"{s}_ReadsPerGene.out.tab"),
        "mapq": os.path.join(st, f"{s}_mapq.txt"),
        "dedup": os.path.join(st, f"{s}_dedupMetrics.txt"),
        "final": os.path.join(st, f"{s}_final_pairs.txt"),
    }


def shared_paths(outdir: str) -> Dict[str, str]:
    c = os.path.join(outdir, "counts")
    return {
        "forward": os.path.join(c, "counts_forward.txt"),
        "reverse": os.path.join(c, "counts_reverse.txt"),
        "deseq2": os.path.join(c, "counts_deseq2.txt"),
        "gene_info": os.path.join(c, "gene_info.tsv"),
    }


# --------------------------------------------------------------------------
# Config / samplesheet
# --------------------------------------------------------------------------

def _merge(a: dict, b: dict) -> dict:
    out = dict(a)
    for k, v in (b or {}).items():
        out[k] = _merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def load_config(paths) -> dict:
    """Merge YAML config files left to right (later files override, dicts merge recursively)."""
    if isinstance(paths, str):
        paths = [paths]
    cfg: dict = {}
    for p in paths or []:
        if p and os.path.exists(p):
            try:
                with open(p) as fh:
                    cfg = _merge(cfg, yaml.safe_load(fh) or {})
            except Exception as exc:
                print(f"qc_summary: could not read config {p}: {exc}", file=sys.stderr)
    return cfg


def read_samplesheet(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = {"Name", "Condition", "Replicate"} - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing column(s): {sorted(missing)}")
    return df


# --------------------------------------------------------------------------
# Parsers
# --------------------------------------------------------------------------

def _num(text: str) -> Optional[float]:
    """Parse '12,345', '85.12%', '0,15' (decimal comma) -> float; None if not numeric."""
    t = text.strip().rstrip("%").strip()
    if not t:
        return None
    if "," in t and "." not in t and re.fullmatch(r"-?\d+,\d+", t) and not re.fullmatch(r"-?\d{1,3}(,\d{3})+", t):
        t = t.replace(",", ".")  # decimal comma (Picard writes the JVM locale)
    else:
        t = t.replace(",", "")
    try:
        v = float(t)
    except ValueError:
        return None
    return None if math.isnan(v) else v


def parse_star_log(path: str) -> Dict[str, Optional[float]]:
    """STAR Log.final.out -> input_pairs, unique_pct, multi_pct, too_short_pct, other_unmapped_pct.

    Lines look like '   Uniquely mapped reads % |\t85.12%'. Missing fields -> None. A file with none of the
    expected labels returns all None.
    """
    vals: Dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            if "|" not in line:
                continue
            key, _, val = line.partition("|")
            v = _num(val)
            if v is not None:
                vals[re.sub(r"\s+", " ", key.strip()).lower()] = v

    def g(label: str) -> Optional[float]:
        return vals.get(label.lower())

    def add(*labels: str) -> Optional[float]:
        got = [g(x) for x in labels]
        got = [x for x in got if x is not None]
        return sum(got) if got else None

    return {
        "input_pairs": g("Number of input reads"),
        "unique_pct": g("Uniquely mapped reads %"),
        "multi_pct": add("% of reads mapped to multiple loci", "% of reads mapped to too many loci"),
        "too_short_pct": g("% of reads unmapped: too short"),
        "other_unmapped_pct": add("% of reads unmapped: too many mismatches", "% of reads unmapped: other"),
    }


def parse_pair_count(path: str) -> Optional[float]:
    """'<sample>\\t<count>' file -> count (last numeric field of the first line that has one)."""
    with open(path) as fh:
        for line in fh:
            for tok in reversed(re.split(r"\s+", line.strip())):
                v = _num(tok)
                if v is not None:
                    return v
    return None


def parse_picard_metrics(path: str) -> Dict[str, Optional[float]]:
    """Picard MarkDuplicates metrics file -> duplicate_pct (PERCENT_DUPLICATION * 100).

    The file has '#' comment lines, a '## METRICS CLASS' line, a tab-separated header line, one data row per
    library, a blank line and then a '## HISTOGRAM' block. Columns are looked up by name. With several
    libraries the value is averaged weighted by READ_PAIRS_EXAMINED. If PERCENT_DUPLICATION is absent it is
    derived from READ_PAIR_DUPLICATES / READ_PAIRS_EXAMINED (documented fallback). Empty/garbled -> None.
    """
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    header: Optional[List[str]] = None
    rows: List[List[str]] = []
    in_metrics = False
    for ln in lines:
        if ln.startswith("## METRICS CLASS"):
            in_metrics, header = True, None
            continue
        if ln.startswith("## HISTOGRAM"):
            break
        if not in_metrics:
            continue
        if not ln.strip():
            if header is not None:
                break
            continue
        if ln.startswith("#"):
            continue
        parts = ln.split("\t")
        if header is None:
            header = parts
        else:
            rows.append(parts)
    out: Dict[str, Optional[float]] = {"duplicate_pct": None}
    if header is None or not rows:
        return out
    col = {name: i for i, name in enumerate(header)}

    def cell(r: List[str], name: str) -> Optional[float]:
        i = col.get(name)
        return _num(r[i]) if i is not None and i < len(r) else None

    pcts, weights = [], []
    for r in rows:
        p = cell(r, "PERCENT_DUPLICATION")
        ex = cell(r, "READ_PAIRS_EXAMINED")
        if p is None:
            d = cell(r, "READ_PAIR_DUPLICATES")
            if d is not None and ex:
                p = d / ex
        if p is not None:
            pcts.append(p * 100.0)
            weights.append(ex if ex else 1.0)
    if pcts:
        out["duplicate_pct"] = float(np.average(pcts, weights=weights)) if sum(weights) > 0 else float(np.mean(pcts))
    return out


def parse_star_genecounts(path: str, strandedness: str) -> Dict[str, Optional[float]]:
    """STAR ReadsPerGene.out.tab -> sense_pct_star.

    Columns: gene_id, unstranded, read 1 same strand as gene, read 1 opposite strand. Rows starting with
    'N_' are STAR summary rows and are skipped. 'forward' => column 3 is sense; 'dUTP' => column 4 is sense.
    """
    same = opp = 0.0
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or p[0].startswith("N_"):
                continue
            a, b = _num(p[2]), _num(p[3])
            if a is None or b is None:
                continue
            same += a
            opp += b
    tot = same + opp
    if tot <= 0:
        return {"sense_pct_star": None}
    sense = opp if strandedness == "dUTP" else same
    return {"sense_pct_star": 100.0 * sense / tot}


def read_matrix(path: str) -> pd.DataFrame:
    """Gene x sample count matrix (tab or comma separated, first column = gene id) -> float DataFrame.

    Non-numeric columns are dropped; duplicated gene ids are summed.
    """
    with open(path) as fh:
        head = fh.readline()
    sep = "\t" if "\t" in head else ","
    df = pd.read_csv(path, sep=sep, index_col=0)
    df.index = df.index.astype(str)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(axis=1, how="all")
    if df.index.duplicated().any():
        df = df.groupby(level=0).sum()
    return df.fillna(0.0)


def read_gene_info(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in df.columns:
        raise ValueError(f"{path}: no gene_id column")
    return df.set_index("gene_id")


# --------------------------------------------------------------------------
# Count-matrix metrics
# --------------------------------------------------------------------------

def log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """log2(CPM + 1) per column; columns with zero total are all zero."""
    lib = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.div(lib, axis=1) * 1e6 + 1.0).fillna(0.0)


def _strip_version(i: str) -> str:
    return re.sub(r"\.\d+$", "", i)


def match_gene_info(gene_ids: pd.Index, info: pd.DataFrame) -> pd.DataFrame:
    """gene_info rows aligned to gene_ids (exact id, else id without '.version'); NaN where unmatched."""
    exact = info.reindex(gene_ids)
    if info.empty or gene_ids.isin(info.index).any():
        return exact
    stripped = info.copy()
    stripped.index = [_strip_version(i) for i in stripped.index]
    stripped = stripped[~stripped.index.duplicated()]
    return stripped.reindex([_strip_version(i) for i in gene_ids]).set_axis(gene_ids)


def replicate_correlations(lcpm: pd.DataFrame, condition_of: Dict[str, str]):
    """Pairwise Pearson r of log2(CPM+1) within each condition.

    Returns (per_sample, pairs): per_sample maps sample -> mean r vs its replicates (absent for singletons),
    pairs is a list of (condition, sample_a, sample_b, r). Genes with zero counts in both samples are excluded
    from that pair. Constant vectors give no value (never 0 or 1).
    """
    pairs = []
    by_cond: Dict[str, List[str]] = {}
    for s in lcpm.columns:
        if s in condition_of:
            by_cond.setdefault(condition_of[s], []).append(s)
    for cond, ss in by_cond.items():
        for i in range(len(ss)):
            for j in range(i + 1, len(ss)):
                a, b = lcpm[ss[i]].to_numpy(), lcpm[ss[j]].to_numpy()
                keep = (a > 0) | (b > 0)
                if keep.sum() < 3 or np.std(a[keep]) == 0 or np.std(b[keep]) == 0:
                    continue
                pairs.append((cond, ss[i], ss[j], float(np.corrcoef(a[keep], b[keep])[0, 1])))
    per: Dict[str, List[float]] = {}
    for _c, a, b, r in pairs:
        per.setdefault(a, []).append(r)
        per.setdefault(b, []).append(r)
    return {s: float(np.mean(v)) for s, v in per.items()}, pairs


def sense_pct_from_strand_matrices(fwd: pd.DataFrame, rev: pd.DataFrame, info: pd.DataFrame) -> Dict[str, float]:
    """Per sample % of counts (on stranded genes) in the strand matching the gene annotation."""
    genes = fwd.index.intersection(rev.index)
    ann = match_gene_info(genes, info)
    strand = ann["strand"] if "strand" in ann.columns else pd.Series(index=genes, dtype=object)
    plus, minus = (strand == "+").to_numpy(), (strand == "-").to_numpy()
    if not (plus.any() or minus.any()):
        return {}
    f, r = fwd.loc[genes], rev.loc[genes]
    out = {}
    for s in f.columns.intersection(r.columns):
        sense = f[s].to_numpy()[plus].sum() + r[s].to_numpy()[minus].sum()
        anti = r[s].to_numpy()[plus].sum() + f[s].to_numpy()[minus].sum()
        if sense + anti > 0:
            out[s] = 100.0 * sense / (sense + anti)
    return out


def rrna_pct_from_matrix(counts: pd.DataFrame, info: pd.DataFrame) -> Dict[str, float]:
    """Per sample % of all counted reads falling in genes with gene_type rRNA / Mt_rRNA."""
    if "gene_type" not in info.columns:
        return {}
    ann = match_gene_info(counts.index, info)
    if ann["gene_type"].isna().all():
        return {}
    is_r = ann["gene_type"].isin(RRNA_TYPES).to_numpy()
    tot = counts.sum(axis=0)
    r = counts[is_r].sum(axis=0)
    return {s: 100.0 * r[s] / tot[s] for s in counts.columns if tot[s] > 0}


# --------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------

def _try(fn, path, label, notes):
    if not os.path.exists(path):
        return None
    try:
        return fn(path)
    except Exception as exc:
        notes.append(f"{label}: could not parse {path}: {exc}")
        return None


def build_sample_qc(outdir: str, sheet: pd.DataFrame, config: dict, counts_path: Optional[str] = None,
                    notes: Optional[List[str]] = None) -> pd.DataFrame:
    """One row per sample with METRIC_COLUMNS; missing inputs give None (blank)."""
    notes = notes if notes is not None else []
    strandedness = str(config.get("strandedness", "dUTP"))
    shared = shared_paths(outdir)
    condition_of = dict(zip(sheet["Name"], sheet["Condition"]))
    recs = {s: {c: None for c in METRIC_COLUMNS} for s in sheet["Name"]}
    for s, cond, rep in zip(sheet["Name"], sheet["Condition"], sheet["Replicate"]):
        r = recs[s]
        r.update(sample=s, condition=cond, replicate=rep)
        p = sample_paths(outdir, s)
        r.update(_try(parse_star_log, p["star_log"], f"{s} STAR log", notes) or {})
        if r["input_pairs"] is not None:
            r["input_pairs_million"] = r["input_pairs"] / 1e6
        r["mapq_pairs"] = _try(parse_pair_count, p["mapq"], f"{s} mapq stats", notes)
        r["final_pairs"] = _try(parse_pair_count, p["final"], f"{s} final stats", notes)
        if r["final_pairs"] is not None and r["input_pairs"]:
            r["final_pct_of_input"] = 100.0 * r["final_pairs"] / r["input_pairs"]
        r.update(_try(parse_picard_metrics, p["dedup"], f"{s} Picard metrics", notes) or {})
        r.update(_try(lambda x: parse_star_genecounts(x, strandedness), p["star_genecnt"], f"{s} ReadsPerGene", notes) or {})

    info = _try(read_gene_info, shared["gene_info"], "gene_info", notes)
    fwd = _try(read_matrix, shared["forward"], "forward counts", notes)
    rev = _try(read_matrix, shared["reverse"], "reverse counts", notes)
    mpath = counts_path if counts_path and os.path.exists(counts_path) else shared["deseq2"]
    mat = _try(read_matrix, mpath, "count matrix", notes)

    if info is not None and fwd is not None and rev is not None:
        for s, v in sense_pct_from_strand_matrices(fwd, rev, info).items():
            if s in recs:
                recs[s]["sense_pct_htseq"] = v
    # total-read matrix for rRNA: forward + reverse when both exist (all counted reads), else the given matrix
    total = None
    if fwd is not None and rev is not None:
        cols = fwd.columns.intersection(rev.columns)
        genes = fwd.index.intersection(rev.index)
        total = fwd.loc[genes, cols] + rev.loc[genes, cols]
    elif mat is not None:
        total = mat
    if total is not None:
        total = total[[c for c in total.columns if c in recs]]
        for s in total.columns:
            recs[s]["counted_reads"] = float(total[s].sum())
        if info is not None:
            for s, v in rrna_pct_from_matrix(total, info).items():
                recs[s]["rrna_pct"] = v
    if mat is not None:
        m = mat[[c for c in mat.columns if c in recs]]
        if m.shape[1] >= 2:
            per, _pairs = replicate_correlations(log_cpm(m), condition_of)
            for s, v in per.items():
                recs[s]["replicate_correlation"] = v
        elif mat.shape[1] and not m.shape[1]:
            notes.append(f"count matrix {mpath}: none of its columns match the samplesheet Names")
    for r in recs.values():
        if r["sense_pct_htseq"] is not None:
            r["sense_pct"], r["sense_source"] = r["sense_pct_htseq"], "htseq"
        elif r["sense_pct_star"] is not None:
            r["sense_pct"], r["sense_source"] = r["sense_pct_star"], "star"
    df = pd.DataFrame([recs[s] for s in sheet["Name"]], columns=METRIC_COLUMNS)
    return df


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Per-sample QC table for TTpipe.")
    ap.add_argument("--outdir", required=True, help="pipeline output directory (outputDir)")
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--config", nargs="*", default=[], help="config YAML file(s); later ones override earlier ones")
    ap.add_argument("--counts", default=None, help="count matrix for replicate correlation (analysis_counts or counts_deseq2.txt)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    notes: List[str] = []
    sheet = read_samplesheet(a.samplesheet)
    df = build_sample_qc(a.outdir, sheet, load_config(a.config), a.counts, notes)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    df.round(4).to_csv(a.out, sep="\t", index=False, na_rep="")
    for n in notes:
        print(f"qc_summary: warning: {n}", file=sys.stderr)
    blank = [c for c in METRIC_COLUMNS[3:] if c != "sense_source" and df[c].isna().all()]
    print(f"qc_summary: wrote {a.out} ({len(df)} sample(s)); no data for: {', '.join(blank) or 'nothing missing'}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
