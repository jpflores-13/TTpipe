#!/usr/bin/env python3
"""Stranded signal plots from the TTpipe bigWigs: one self-contained interactive HTML page.

    python scripts/stranded_signal.py \
        --bigwig DMSO=forward=OUT/tracks/DMSO_forward.bw --bigwig DMSO=reverse=OUT/tracks/DMSO_reverse.bw ... \
        --condition-order DMSO,2H --scale-factors OUT/tracks/scale_factors.tsv \
        --gene-info OUT/counts/gene_info.tsv --gtf genes.gtf \
        --counts-forward OUT/counts/counts_forward.txt --counts-reverse OUT/counts/counts_reverse.txt \
        --counts-matrix OUT/counts/counts_deseq2.txt --samples samplesheet.txt \
        [--loci Gapdh chr1:1000-5000 ...] [--degs OUT/analysis/deseq2/degs_*.csv --top-degs 10] [--regions-bed r.bed] \
        --out OUT/analysis/TTpipe_stranded_signal.html

The page has three parts (plotly.js is inlined, nothing is fetched, it works from file://):
  1. Locus tracks: per selected locus and condition, the + strand RNA above the axis and the - strand RNA mirrored
     below it, a shared x axis, a y scale shared across conditions (or free per condition), and gene models.
  2. Stranded metagene profiles: sense and antisense signal relative to gene orientation (TSS-centred, and gene
     body scaled to a fixed number of bins), per gene_type group and condition.
  3. A strand table: the antisense fraction of the read counts per gene and condition.

The module separates data preparation (functions returning plain dicts / arrays, tested without a browser) from
the HTML templating (build_html). The bigWigs are read with pyBigWig directly (not deepTools).

Conventions
  forward = RNA transcribed from the + genomic strand, reverse = from the - strand. Both bigWigs hold POSITIVE values
  (bamCoverage), so the reverse strand is mirrored here. A + gene's sense signal is in `forward`, a - gene's in
  `reverse`. gene_info start/end are the GTF's 1-based inclusive coordinates (scripts/gene_info.py); bigWig
  coordinates are 0-based half-open, so a gene covers [start - 1, end) in bigWig coordinates.
  The scale factors are already applied inside the bigWigs: nothing is rescaled here.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import gzip
import json
import math
import os
import re
import sys
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gene_info import parse_attributes  # noqa: E402
from results_explorer import CSS as BASE_CSS, sig_round, to_json  # noqa: E402

MAX_BINS = 4000            # a locus window never gets more bins than this (the bin size is raised in whole multiples)
MIN_GROUP_GENES = 20       # metagene groups with fewer genes are skipped
# GENCODE merged these gene_types into "lncRNA" (human v31 / mouse vM22 onward, applied unevenly in some
# releases, e.g. vM25 still says lincRNA / antisense). A metagene group named "lncRNA" also takes them, so
# older annotations do not come out with an empty lncRNA group.
LNCRNA_LEGACY_TYPES = frozenset({"lncrna", "lincrna", "antisense", "antisense_rna", "sense_intronic",
                                 "sense_overlapping", "processed_transcript", "3prime_overlapping_ncrna",
                                 "bidirectional_promoter_lncrna", "macro_lncrna", "non_coding"})
TSS_BIN_BP = 50            # bin size of the TSS-centred metagene
FLANK_BIN_BP = 100         # bin size of the flanks in the gene-body metagene
BODY_BINS = 100            # bins across the scaled gene body
MAX_MODEL_ROWS = 10        # gene-model rows drawn per locus
SIZE_WARN_MB = 40.0
ALLOWED_STATIC = ("svg", "pdf", "png")
NORMALIZATION = ("scaled coverage: the bigWigs were written by bamCoverage with the scale factors of "
                 "tracks/scale_factors.tsv already applied (nothing is rescaled here); values are averaged in bins. "
                 "The factors are library-size factors, or spike-in factors when spike_in.normalize includes "
                 "'tracks' (the table's source column says which).")


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------

def _open_text(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def log(msg: str) -> None:
    print(f"stranded_signal: {msg}", file=sys.stderr)


def _import_pybigwig():
    try:
        import pyBigWig  # noqa: N813
    except ImportError as e:  # pragma: no cover - environment dependent
        raise SystemExit("stranded_signal needs the Python package pyBigWig (conda: bioconda::pybigwig, pip: pybigwig); "
                         f"it is not importable here: {e}")
    return pyBigWig


_VER_RE = re.compile(r"^([A-Za-z]+\d+)\.\d+$")


def _strip_version(gene_id: str) -> str:
    """ENSG00000123456.12 -> ENSG00000123456 (only Ensembl/GENCODE-looking ids; 'Gene.1' style names are kept)."""
    m = _VER_RE.match(gene_id)
    return m.group(1) if m else gene_id


# --------------------------------------------------------------------------
# Input tables
# --------------------------------------------------------------------------

def read_gene_table(path: str) -> pd.DataFrame:
    """gene_info.tsv -> DataFrame with gene_id, gene_name, gene_type, chr, start (1-based), end, strand, length."""
    gi = pd.read_csv(path, sep="\t", dtype={"gene_id": str, "gene_name": str, "gene_type": str, "chr": str,
                                              "strand": str}, keep_default_na=False)
    need = {"gene_id", "chr", "start", "end", "strand"}
    if not need <= set(gi.columns):
        raise ValueError(f"{path}: gene table needs columns {sorted(need)}; has {list(gi.columns)}")
    for c in ("gene_name", "gene_type"):
        if c not in gi.columns:
            gi[c] = ""
    gi["gene_name"] = np.where(gi["gene_name"].str.strip() == "", gi["gene_id"], gi["gene_name"])
    gi["start"] = pd.to_numeric(gi["start"]).astype(int)
    gi["end"] = pd.to_numeric(gi["end"]).astype(int)
    return gi.drop_duplicates("gene_id").reset_index(drop=True)


def read_scale_factors(path: Optional[str]) -> Dict[str, float]:
    if not path or not os.path.exists(path):
        return {}
    sf = pd.read_csv(path, sep="\t")
    if not {"condition", "scale_factor"} <= set(sf.columns):
        return {}
    return {str(r.condition): float(r.scale_factor) for r in sf.itertuples()}


def read_bed(path: str) -> List[Dict[str, Any]]:
    """BED3-BED6 -> list of dicts (chrom, start0, end0, name or '', strand '+', '-' or '')."""
    out = []
    with _open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t") if "\t" in line else line.split()
            if len(f) < 3:
                continue
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                continue
            if e <= s:
                continue
            name = f[3] if len(f) > 3 and f[3] not in (".", "") else ""
            strand = f[5] if len(f) > 5 and f[5] in ("+", "-") else ""
            out.append({"chrom": f[0], "start0": s, "end0": e, "name": name, "strand": strand})
    return out


def contrast_name_from_path(path: str) -> str:
    m = re.fullmatch(r"degs_(.+)\.csv", os.path.basename(path))
    return m.group(1) if m else os.path.splitext(os.path.basename(path))[0]


def top_deg_ids(path: str, n: int) -> List[str]:
    """The n gene_ids with the smallest non-NA padj of one DESeq2 table (ties keep file order)."""
    df = pd.read_csv(path)
    if df.shape[1] == 0:
        return []
    ids = df.iloc[:, 0].astype(str)
    padj = pd.to_numeric(df["padj"], errors="coerce") if "padj" in df.columns else pd.Series(np.nan, index=df.index)
    ok = padj.notna()
    order = padj[ok].sort_values(kind="mergesort").index[:n]
    return [ids[i] for i in order]


# --------------------------------------------------------------------------
# Locus selection
# --------------------------------------------------------------------------

COORD_RE = re.compile(r"^([^\s:]+):([\d,]+)-([\d,]+)$")


class GeneIndex:
    """Lookups over the gene table: by gene_id (exact, then without version) and by symbol (case-insensitive)."""

    def __init__(self, gi: pd.DataFrame):
        self.gi = gi
        self.by_id: Dict[str, int] = {}
        self.by_id_nov: Dict[str, List[int]] = {}
        self.by_name: Dict[str, List[int]] = {}
        for i, (gid, nm) in enumerate(zip(gi["gene_id"], gi["gene_name"])):
            self.by_id[gid] = i
            self.by_id_nov.setdefault(_strip_version(gid), []).append(i)
            self.by_name.setdefault(str(nm).lower(), []).append(i)

    def lookup(self, entry: str) -> Tuple[List[int], str]:
        """(row indices, kind) with kind in 'id', 'id_noversion', 'symbol', 'ambiguous_symbol' or '' (not found)."""
        if entry in self.by_id:
            return [self.by_id[entry]], "id"
        nov = _strip_version(entry)
        if nov in self.by_id_nov:
            return list(self.by_id_nov[nov]), "id_noversion"
        hits = self.by_name.get(entry.lower(), [])
        if len(hits) == 1:
            return list(hits), "symbol"
        if len(hits) > 1:
            return list(hits), "ambiguous_symbol"
        return [], ""


def gene_locus(row: pd.Series, source: str) -> Dict[str, Any]:
    return {"key": ("gene", row["gene_id"]), "label": row["gene_name"], "chrom": row["chr"],
            "start0": int(row["start"]) - 1, "end0": int(row["end"]), "strand": row["strand"],
            "gene_id": row["gene_id"], "gene_name": row["gene_name"], "sources": [source]}


def region_locus(chrom: str, start0: int, end0: int, label: str, source: str, strand: str = "") -> Dict[str, Any]:
    return {"key": ("region", chrom, start0, end0), "label": label or f"{chrom}:{start0 + 1}-{end0}",
            "chrom": chrom, "start0": start0, "end0": end0, "strand": strand, "gene_id": "", "gene_name": "",
            "sources": [source]}


def select_loci(entries: Sequence[str], gene_table: pd.DataFrame, deg_tables: Optional[Dict[str, str]] = None,
                top_n: int = 0, bed: Optional[str] = None, max_loci: int = 40
                ) -> Tuple[List[Dict[str, Any]], List[Dict[str, str]], List[str]]:
    """Resolve the configured loci. Returns (loci, not_found, notes).

    Order: `track_loci` entries, then the top DEGs of each contrast, then the BED regions. Entries are combined,
    de-duplicated (a gene reached twice is listed once with all its sources) and capped at `max_loci`. Entries that
    cannot be resolved go to `not_found` (never an exception)."""
    idx = GeneIndex(gene_table)
    loci: List[Dict[str, Any]] = []
    seen: Dict[tuple, Dict[str, Any]] = {}
    not_found: List[Dict[str, str]] = []
    notes: List[str] = []

    def add(loc: Dict[str, Any]) -> None:
        if loc["key"] in seen:
            for s in loc["sources"]:
                if s not in seen[loc["key"]]["sources"]:
                    seen[loc["key"]]["sources"].append(s)
            return
        seen[loc["key"]] = loc
        loci.append(loc)

    for entry in entries or []:
        entry = str(entry).strip()
        if not entry:
            continue
        m = COORD_RE.match(entry)
        if m:
            s1, e1 = int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
            if s1 < 1 or e1 < s1:
                not_found.append({"entry": entry, "reason": "coordinates must be chr:start-end with 1 <= start <= end"})
                continue
            add(region_locus(m.group(1), s1 - 1, e1, f"{m.group(1)}:{s1}-{e1}", "coordinates"))
            continue
        rows, kind = idx.lookup(entry)
        if not rows:
            not_found.append({"entry": entry, "reason": "no gene with this symbol or gene_id in gene_info.tsv"})
            continue
        if kind == "ambiguous_symbol":
            ids = ", ".join(gene_table.at[i, "gene_id"] for i in rows)
            notes.append(f"symbol {entry!r} matches {len(rows)} genes ({ids}): all are shown")
        for i in rows:
            loc = gene_locus(gene_table.iloc[i], "config: " + entry)
            if kind == "ambiguous_symbol":
                loc["label"] = f"{loc['gene_name']} ({loc['gene_id']})"
            add(loc)

    if top_n and top_n > 0:
        for cname, path in (deg_tables or {}).items():
            if not path or not os.path.exists(path):
                continue
            for gid in top_deg_ids(path, top_n):
                rows, kind = idx.lookup(gid)
                if not rows:
                    not_found.append({"entry": f"{gid} (top DEG, {cname})", "reason": "gene_id not in gene_info.tsv"})
                    continue
                for i in rows[:1]:
                    add(gene_locus(gene_table.iloc[i], f"top DEG: {cname}"))

    if bed:
        if not os.path.exists(bed):
            not_found.append({"entry": bed, "reason": "track_regions_bed file does not exist"})
        else:
            regs = read_bed(bed)
            if not regs:
                not_found.append({"entry": bed, "reason": "no usable BED lines (need chrom, start, end)"})
            for r in regs:
                add(region_locus(r["chrom"], r["start0"], r["end0"], r["name"], "BED regions", r["strand"]))

    if len(loci) > max_loci:
        notes.append(f"{len(loci)} loci were selected; only the first {max_loci} are drawn (track_max_loci). "
                     "Order: track_loci, top DEGs per contrast, BED regions.")
        loci = loci[:max_loci]
    return loci, not_found, notes


# --------------------------------------------------------------------------
# bigWig binning
# --------------------------------------------------------------------------

def plan_bins(start0: int, end0: int, chrom_len: int, bin_size: int, max_bins: int = MAX_BINS) -> Tuple[int, int, int]:
    """(window start, n bins, bin width) for a window: whole bins of `bin_size` (raised in whole multiples when the
    window would need more than `max_bins`), never past the chromosome end."""
    length = end0 - start0
    bs = int(bin_size)
    if length > bs * max_bins:
        bs = int(bin_size) * math.ceil(length / (int(bin_size) * max_bins))
    n = max(1, math.ceil(length / bs))
    if start0 + n * bs > chrom_len:
        n = (chrom_len - start0) // bs
        if n < 1:
            n, bs = 1, max(1, chrom_len - start0)
    return start0, n, bs


def bin_mean(bw, chrom: str, start0: int, end0: int, nbins: int) -> np.ndarray:
    """Mean signal per bin over [start0, end0) split into `nbins` equal bins, uncovered bases counted as 0.

    The per-base values are read once (pyBigWig `values`, NaN where the bigWig has no entry, NaN -> 0) and summed per
    bin with a cumulative sum; bin edges are rounded to whole bases (a gene body of 2,001 bp in 100 bins has bins of
    20 or 21 bp). This is about 100 times faster than pyBigWig's `stats(nBins=...)`, which also differs in meaning:
    its `mean` ignores uncovered bases (a half-covered bin returns the value of the covered half)."""
    start0, end0 = int(start0), int(end0)
    try:
        v = bw.values(chrom, start0, end0, numpy=True)
    except TypeError:                                  # pyBigWig built without numpy support
        v = bw.values(chrom, start0, end0)
    v = np.nan_to_num(np.asarray(v, dtype=float), nan=0.0)
    edges = np.rint(np.linspace(0, v.size, int(nbins) + 1)).astype(np.int64)
    cs = np.concatenate(([0.0], np.cumsum(v)))
    widths = np.maximum(np.diff(edges), 1)
    return (cs[edges[1:]] - cs[edges[:-1]]) / widths


def locus_signal(bws: Dict[str, Dict[str, Any]], conditions: Sequence[str], chrom: str, start0: int, end0: int,
                 bin_size: int) -> Optional[Dict[str, Any]]:
    """Binned forward / reverse signal for one window and all conditions. None if the chromosome is not in the
    bigWigs. Window is clipped to the chromosome."""
    first = bws[conditions[0]]["forward"]
    clen = first.chroms(chrom)
    if not clen:
        return None
    start0 = max(0, start0)
    end0 = min(int(clen), end0)
    if end0 <= start0:
        return None
    s, n, bs = plan_bins(start0, end0, int(clen), bin_size)
    sig = {}
    for c in conditions:
        sig[c] = {"f": bin_mean(bws[c]["forward"], chrom, s, s + n * bs, n),
                  "r": bin_mean(bws[c]["reverse"], chrom, s, s + n * bs, n)}
    return {"start0": s, "n": n, "bin": bs, "sig": sig}


# --------------------------------------------------------------------------
# Gene models from the GTF
# --------------------------------------------------------------------------

def read_gene_models(gtf: str, windows: Dict[str, List[Tuple[int, int]]]) -> Dict[str, Dict[str, Any]]:
    """Stream the GTF once and keep only exons (and gene lines) that overlap a window.

    `windows`: chrom -> list of (start1, end1), 1-based inclusive. Returns gene_id -> {gene_id, gene_name, chrom,
    strand, exons: [(s1, e1), ...], extent: (s1, e1) or None}. Lines on other chromosomes or outside the windows are
    skipped before their attributes are parsed, so a very large GTF costs little memory."""
    genes: Dict[str, Dict[str, Any]] = {}
    if not windows:
        return genes
    with _open_text(gtf) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.split("\t", 8)
            if len(f) < 9:
                continue
            wins = windows.get(f[0])
            if not wins:
                continue
            feat = f[2]
            if feat != "exon" and feat != "gene":
                continue
            try:
                s, e = int(f[3]), int(f[4])
            except ValueError:
                continue
            if not any(s <= we and e >= ws for ws, we in wins):
                continue
            a = parse_attributes(f[8])
            gid = a.get("gene_id")
            if not gid:
                continue
            g = genes.setdefault(gid, {"gene_id": gid, "gene_name": a.get("gene_name") or gid, "chrom": f[0],
                                       "strand": f[6] if f[6] in ("+", "-") else ".", "exons": [], "extent": None})
            if a.get("gene_name"):
                g["gene_name"] = a["gene_name"]
            if feat == "gene":
                g["extent"] = (s, e)
            else:
                g["exons"].append((s, e))
    return genes


def merge_intervals(iv: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    out: List[List[int]] = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def models_for_window(genes: Dict[str, Dict[str, Any]], gene_extent: Dict[str, Tuple[int, int]], chrom: str,
                      w_start1: int, w_end1: int, max_rows: int = MAX_MODEL_ROWS) -> Tuple[List[Dict[str, Any]], int]:
    """Collapsed gene models (union of all exons of a gene, 1-based inclusive) for genes overlapping the window,
    assigned to non-overlapping display rows. Returns (models, number of genes not drawn for lack of rows)."""
    ms = []
    for g in genes.values():
        if g["chrom"] != chrom:
            continue
        ex = merge_intervals(g["exons"])
        ext = g.get("extent") or gene_extent.get(g["gene_id"])
        if ext is None and ex:
            ext = (ex[0][0], ex[-1][1])
        if ext is None or ext[1] < w_start1 or ext[0] > w_end1:
            continue
        ex = [x for x in ex if x[1] >= w_start1 and x[0] <= w_end1]
        ms.append({"id": g["gene_id"], "name": g["gene_name"], "strand": g["strand"], "start": ext[0], "end": ext[1],
                   "exons": [[s, e] for s, e in ex]})
    ms.sort(key=lambda m: (m["start"], m["end"], m["id"]))
    pad = max(1, (w_end1 - w_start1) // 40)       # keep room for the gene label
    row_end: List[int] = []
    drawn, hidden = [], 0
    for m in ms:
        for r, last in enumerate(row_end):
            if m["start"] > last + pad:
                row_end[r] = m["end"]
                m["row"] = r
                break
        else:
            if len(row_end) >= max_rows:
                hidden += 1
                continue
            row_end.append(m["end"])
            m["row"] = len(row_end) - 1
        drawn.append(m)
    return drawn, hidden


# --------------------------------------------------------------------------
# Locus tracks
# --------------------------------------------------------------------------

def sig_list(a: np.ndarray, digits: int = 3) -> List[Optional[float]]:
    return [sig_round(float(v), digits) for v in a]


def build_locus_tracks(loci: List[Dict[str, Any]], bws: Dict[str, Dict[str, Any]], conditions: Sequence[str],
                       flank: int, bin_size: int, gtf: Optional[str], gene_table: Optional[pd.DataFrame]
                       ) -> Tuple[List[Dict[str, Any]], List[Dict[str, str]]]:
    """Signal and gene models for each locus. Returns (locus records for the page, loci that could not be drawn)."""
    out, failed = [], []
    wins: Dict[str, List[Tuple[int, int]]] = {}
    prepared = []
    for loc in loci:
        sg = locus_signal(bws, conditions, loc["chrom"], loc["start0"] - flank, loc["end0"] + flank, bin_size)
        if sg is None:
            failed.append({"entry": loc["label"], "reason": f"chromosome {loc['chrom']!r} is not in the bigWigs "
                                                           "(or the window is empty)"})
            continue
        prepared.append((loc, sg))
        wins.setdefault(loc["chrom"], []).append((sg["start0"] + 1, sg["start0"] + sg["n"] * sg["bin"]))
    genes = read_gene_models(gtf, wins) if gtf and wins else {}
    extent = {}
    if gene_table is not None:
        extent = {g: (int(s), int(e)) for g, s, e in zip(gene_table["gene_id"], gene_table["start"], gene_table["end"])}
    for loc, sg in prepared:
        w1, w2 = sg["start0"] + 1, sg["start0"] + sg["n"] * sg["bin"]
        models, hidden = models_for_window(genes, extent, loc["chrom"], w1, w2)
        out.append({
            "label": loc["label"], "chrom": loc["chrom"], "start": loc["start0"] + 1, "end": loc["end0"],
            "strand": loc["strand"] or ".", "gene_id": loc["gene_id"], "sources": loc["sources"],
            "win_start": w1, "win_end": w2, "bin": sg["bin"], "n": sg["n"],
            "bin_raised": sg["bin"] != bin_size and sg["bin"] % bin_size == 0,
            "sig": {c: {"f": sig_list(sg["sig"][c]["f"]), "r": sig_list(sg["sig"][c]["r"])} for c in conditions},
            "models": models, "models_hidden": hidden,
        })
    return out, failed


# --------------------------------------------------------------------------
# Metagene profiles
# --------------------------------------------------------------------------

def read_matrix(path: str) -> pd.DataFrame:
    m = pd.read_csv(path, sep="\t", index_col=0)
    m.index = m.index.astype(str)
    m.columns = [str(c) for c in m.columns]
    return m


def expressed_gene_ids(counts_path: Optional[str], min_count: float, min_samples: int) -> Optional[set]:
    """Gene ids with counts >= min_count in >= min_samples samples; None when no matrix is given (all genes)."""
    if not counts_path or not os.path.exists(counts_path):
        return None
    m = read_matrix(counts_path)
    keep = (m.to_numpy(dtype=float) >= float(min_count)).sum(axis=1) >= int(min_samples)
    return set(m.index[keep])


def window_matrix(path: str, chroms: Sequence[str], starts: Sequence[int], ends: Sequence[int], nbins: int) -> np.ndarray:
    """One row per window (genomic left-to-right) of bin means; module-level so process pools can pickle it."""
    pyBigWig = _import_pybigwig()
    bw = pyBigWig.open(path)
    try:
        m = np.zeros((len(chroms), nbins))
        for i, (c, s, e) in enumerate(zip(chroms, starts, ends)):
            m[i] = bin_mean(bw, c, s, e, nbins)
        return m
    finally:
        bw.close()


def valid_windows(chrom_sizes: Dict[str, int], chroms: Sequence[str], starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    return np.array([c in chrom_sizes and s >= 0 and e <= chrom_sizes[c] and e > s
                     for c, s, e in zip(chroms, starts, ends)], dtype=bool)


def _segments(view: str, regions: pd.DataFrame, flank: int) -> Tuple[List[Tuple[np.ndarray, np.ndarray, int]], Dict[str, Any]]:
    """Genomic windows (per region) and bin counts that make one metagene view, plus x-axis information.

    tss  : one window [anchor - flank, anchor + flank) where anchor is the TSS (or the region centre) after
           orientation, TSS_BIN_BP bins.
    body : upstream flank | BODY_BINS bins over the gene | downstream flank (in genomic left-to-right order)."""
    st = regions["start0"].to_numpy()
    en = regions["end0"].to_numpy()
    strand = regions["strand_int"].to_numpy()
    if view == "tss":
        anchor = np.where(strand < 0, en - 1, st) if "anchor" not in regions else regions["anchor"].to_numpy()
        s = np.where(strand < 0, anchor - flank + 1, anchor - flank)
        n = int(round(2 * flank / TSS_BIN_BP))
        x = (-flank + (np.arange(n) + 0.5) * (2 * flank / n)).tolist()
        info = {"x": x, "ticks": {"vals": [-flank, 0, flank], "text": [f"-{flank / 1000:g} kb", "TSS", f"+{flank / 1000:g} kb"]}}
        return [(s, s + 2 * flank, n)], info
    nfl = max(1, int(round(flank / FLANK_BIN_BP)))
    segs = [(st - flank, st, nfl), (st, en, BODY_BINS), (en, en + flank, nfl)]
    total = 2 * nfl + BODY_BINS
    info = {"x": list(range(total)),
            "ticks": {"vals": [0, nfl, nfl + BODY_BINS, total - 1],
                      "text": [f"-{flank / 1000:g} kb", "TSS", "TES", f"+{flank / 1000:g} kb"]}}
    return segs, info


def orient_profiles(fmat: np.ndarray, rmat: np.ndarray, strand_int: np.ndarray, oriented: bool
                    ) -> Tuple[np.ndarray, np.ndarray]:
    """From genomic-orientation forward / reverse matrices to (sense, antisense), 5' on the left.

    + gene: sense = forward, antisense = reverse (columns kept). - gene: sense = reverse, antisense = forward, and the
    columns are flipped so that the 5' end is still on the left. With oriented=False (a region set without strand)
    nothing is swapped or flipped: 'sense' is the + strand RNA and 'antisense' the - strand RNA."""
    if not oriented:
        return fmat, rmat
    neg = (strand_int < 0)[:, None]
    sense = np.where(neg, rmat[:, ::-1], fmat)
    anti = np.where(neg, fmat[:, ::-1], rmat)
    return sense, anti


def _task(args):
    return window_matrix(*args)


def _run_matrix_tasks(tasks: List[tuple], threads: int) -> List[np.ndarray]:
    if threads > 1 and len(tasks) > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=threads) as ex:
            return list(ex.map(_task, tasks))
    return [_task(t) for t in tasks]


def profile_group(bigwigs: Dict[str, Dict[str, str]], conditions: Sequence[str], regions: pd.DataFrame, flank: int,
                  view: str, oriented: bool, threads: int = 1) -> Tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, Any], int]:
    """Per condition sense / antisense matrices (regions x bins) for one view. Regions whose window leaves the
    chromosome (or whose chromosome is missing from the bigWigs) are dropped; returns (matrices, x info, n kept)."""
    pyBigWig = _import_pybigwig()
    first = bigwigs[conditions[0]]["forward"]
    bw = pyBigWig.open(first)
    sizes = dict(bw.chroms())
    bw.close()
    segs, info = _segments(view, regions, flank)
    ok = np.ones(len(regions), dtype=bool)
    chroms = regions["chrom"].tolist()
    for s, e, _ in segs:
        ok &= valid_windows(sizes, chroms, s, e)
    keep = np.flatnonzero(ok)
    if keep.size == 0:
        return {}, info, 0
    ch = [chroms[i] for i in keep]
    tasks, order = [], []
    for c in conditions:
        for strand in ("forward", "reverse"):
            for k, (s, e, nb) in enumerate(segs):
                tasks.append((bigwigs[c][strand], ch, s[keep].tolist(), e[keep].tolist(), nb))
                order.append((c, strand, k))
    res = _run_matrix_tasks(tasks, threads)
    parts: Dict[Tuple[str, str], List[np.ndarray]] = {}
    for (c, strand, k), m in zip(order, res):
        parts.setdefault((c, strand), []).append(m)
    sint = regions["strand_int"].to_numpy()[keep]
    out = {}
    for c in conditions:
        fm = np.hstack(parts[(c, "forward")])
        rm = np.hstack(parts[(c, "reverse")])
        se, an = orient_profiles(fm, rm, sint, oriented)
        out[c] = {"sense": se, "anti": an}
    return out, info, int(keep.size)


def summarize_profiles(mats: Dict[str, Dict[str, np.ndarray]]) -> Dict[str, Dict[str, List[Optional[float]]]]:
    return {c: {"sense_mean": sig_list(m["sense"].mean(axis=0), 4), "anti_mean": sig_list(m["anti"].mean(axis=0), 4),
                "sense_median": sig_list(np.median(m["sense"], axis=0), 4),
                "anti_median": sig_list(np.median(m["anti"], axis=0), 4)} for c, m in mats.items()}


def gene_regions(gene_table: pd.DataFrame) -> pd.DataFrame:
    r = pd.DataFrame({"gene_id": gene_table["gene_id"], "chrom": gene_table["chr"],
                      "start0": gene_table["start"] - 1, "end0": gene_table["end"],
                      "strand_int": np.where(gene_table["strand"] == "-", -1, 1)})
    r.loc[~gene_table["strand"].isin(["+", "-"]).to_numpy(), "strand_int"] = 1
    return r.reset_index(drop=True)


def build_metagene(bigwigs: Dict[str, Dict[str, str]], conditions: Sequence[str], gene_table: pd.DataFrame,
                   groups: Sequence[str], flank: int, min_length: int, expressed: Optional[set],
                   bed_regions: Optional[List[Dict[str, Any]]] = None, threads: int = 1
                   ) -> Dict[str, Any]:
    """Stranded metagene profiles for each gene_type group (and a BED group). See the module docstring."""
    min_length = max(int(min_length), BODY_BINS)
    out_groups: List[Dict[str, Any]] = []
    notes: List[str] = []
    gt = gene_table.copy()
    gt["gene_type"] = gt["gene_type"].astype(str)
    # genes on a strand '.' have no orientation: they are left out (a GTF gene line always has + or -)
    gt = gt[gt["strand"].isin(["+", "-"])]
    if expressed is not None:
        gt = gt[gt["gene_id"].isin(expressed)]
    for name in groups:
        key = str(name).lower()
        label = name
        if key == "all":
            sel = gt
        elif key == "lncrna":
            sel = gt[gt["gene_type"].str.lower().isin(LNCRNA_LEGACY_TYPES)]
            legacy = sorted(set(sel["gene_type"]) - {t for t in sel["gene_type"] if t.lower() == "lncrna"})
            if legacy:
                label = f"{name} (incl. {', '.join(legacy)})"
                notes.append(f"group {name!r}: older GENCODE gene_types counted as lncRNA: {', '.join(legacy)}")
        else:
            sel = gt[gt["gene_type"].str.lower() == key]
        n_all = len(sel)
        if n_all < MIN_GROUP_GENES:
            notes.append(f"group {name!r}: only {n_all} expressed gene(s) (< {MIN_GROUP_GENES}), skipped")
            continue
        reg = gene_regions(sel)
        t0 = time.time()
        tm, tinfo, n_tss = profile_group(bigwigs, conditions, reg, flank, "tss", True, threads)
        grp: Dict[str, Any] = {"name": name, "label": label, "kind": "gene", "oriented": True, "n_tss": n_tss,
                               "n_expressed": n_all, "tss": None, "body": None, "body_note": ""}
        if n_tss >= MIN_GROUP_GENES:
            grp["tss"] = {"x": tinfo["x"], "ticks": tinfo["ticks"], "profiles": summarize_profiles(tm)}
        else:
            notes.append(f"group {name!r}: only {n_tss} gene(s) have a full TSS window inside the chromosomes; no TSS view")
        long_reg = reg[(reg["end0"] - reg["start0"]) >= min_length].reset_index(drop=True)
        bm, binfo, n_body = ({}, {}, 0)
        if len(long_reg) >= MIN_GROUP_GENES:
            bm, binfo, n_body = profile_group(bigwigs, conditions, long_reg, flank, "body", True, threads)
        grp["n_body"] = n_body
        if n_body >= MIN_GROUP_GENES:
            grp["body"] = {"x": binfo["x"], "ticks": binfo["ticks"], "profiles": summarize_profiles(bm)}
        else:
            grp["body_note"] = (f"fewer than {MIN_GROUP_GENES} expressed genes of at least {min_length} bp: "
                                "no gene-body view for this group")
            notes.append(f"group {name!r}: {grp['body_note']}")
        log(f"metagene group {name!r}: {n_all} expressed genes; TSS view {n_tss}, body view {n_body} genes "
            f"({time.time() - t0:.1f} s)")
        if grp["tss"] or grp["body"]:
            out_groups.append(grp)
    if bed_regions:
        df = pd.DataFrame(bed_regions)
        oriented = bool((df["strand"].isin(["+", "-"])).all())
        df["strand_int"] = np.where((df["strand"] == "-") & oriented, -1, 1)
        df["anchor"] = (df["start0"] + df["end0"]) // 2
        if len(df) < MIN_GROUP_GENES:
            notes.append(f"metagene_bed: only {len(df)} region(s) (< {MIN_GROUP_GENES}), skipped")
        else:
            t0 = time.time()
            reg = df[["chrom", "start0", "end0", "strand_int", "anchor"]].copy()
            tm, tinfo, n = profile_group(bigwigs, conditions, reg, flank, "tss", oriented, threads)
            if n >= MIN_GROUP_GENES:
                tinfo["ticks"]["text"][1] = "centre"
                out_groups.append({"name": "BED regions", "label": "BED regions (centre-aligned)", "kind": "bed",
                                   "oriented": oriented, "n_tss": n, "n_body": 0, "n_expressed": len(df),
                                   "tss": {"x": tinfo["x"], "ticks": tinfo["ticks"], "profiles": summarize_profiles(tm)},
                                   "body": None, "body_note": "region sets have no gene-body view"})
                if not oriented:
                    notes.append("metagene_bed has no strand column (or not for every region): both strands are shown "
                                 "as + strand RNA (solid) and - strand RNA (dashed), not as sense / antisense")
                log(f"metagene BED group: {n} of {len(df)} regions ({time.time() - t0:.1f} s)")
            else:
                notes.append(f"metagene_bed: only {n} region(s) have a full window inside the chromosomes, skipped")
    return {"groups": out_groups, "notes": notes, "flank": int(flank), "min_length": int(min_length),
            "tss_bin": TSS_BIN_BP, "body_bins": BODY_BINS}


# --------------------------------------------------------------------------
# Strand table
# --------------------------------------------------------------------------

def strand_table(counts_f: pd.DataFrame, counts_r: pd.DataFrame, gene_table: pd.DataFrame,
                 sample_condition: Dict[str, str], conditions: Sequence[str], min_reads: int, top: int = 5000
                 ) -> Tuple[pd.DataFrame, int]:
    """Antisense fraction per gene and condition = antisense / (sense + antisense), counts summed over the samples
    of the condition. sense = forward counts for a + gene, reverse counts for a - gene (antisense the other strand).

    Returns (table, n_genes_passing_min_reads). Genes need >= min_reads total reads over all samples and both
    strands; only the `top` genes by total reads are kept. A condition with sense + antisense = 0 gets NaN."""
    ids = [g for g in gene_table["gene_id"] if g in counts_f.index and g in counts_r.index]
    gi = gene_table.set_index("gene_id").loc[ids]
    f = counts_f.loc[ids]
    r = counts_r.loc[ids]
    strand_plus = (gi["strand"] != "-").to_numpy()[:, None]
    samples = [s for s in f.columns if s in r.columns]
    sense_s = np.where(strand_plus, f[samples].to_numpy(float), r[samples].to_numpy(float))
    anti_s = np.where(strand_plus, r[samples].to_numpy(float), f[samples].to_numpy(float))
    total = sense_s.sum(axis=1) + anti_s.sum(axis=1)
    df = pd.DataFrame({"gene_id": ids, "gene_name": gi["gene_name"].to_numpy(), "gene_type": gi["gene_type"].to_numpy(),
                       "chr": gi["chr"].to_numpy(), "start": gi["start"].to_numpy(), "end": gi["end"].to_numpy(),
                       "strand": gi["strand"].to_numpy(), "total_reads": total})
    for c in conditions:
        cols = [i for i, s in enumerate(samples) if sample_condition.get(s) == c]
        se = sense_s[:, cols].sum(axis=1) if cols else np.zeros(len(ids))
        an = anti_s[:, cols].sum(axis=1) if cols else np.zeros(len(ids))
        den = se + an
        df[f"{c}__sense"] = se
        df[f"{c}__antisense"] = an
        df[f"{c}__antisense_fraction"] = np.where(den > 0, an / np.where(den > 0, den, 1), np.nan)
    df = df[df["total_reads"] >= min_reads]
    n_pass = len(df)
    df = df.sort_values(["total_reads", "gene_id"], ascending=[False, True], kind="mergesort").head(top)
    return df.reset_index(drop=True), n_pass


def strand_table_block(df: pd.DataFrame, conditions: Sequence[str], n_pass: int) -> Dict[str, Any]:
    """Compact rows for the page: [gene_id, name, type, chr, start, end, strand, total, [[sense, anti, fraction] per condition]]."""
    rows = []
    for i in range(len(df)):
        r = df.iloc[i]
        cond = [[int(round(r[f"{c}__sense"])), int(round(r[f"{c}__antisense"])), sig_round(r[f"{c}__antisense_fraction"], 4)]
                for c in conditions]
        rows.append([r["gene_id"], r["gene_name"], r["gene_type"], r["chr"], int(r["start"]), int(r["end"]), r["strand"],
                     int(round(r["total_reads"])), cond])
    return {"conditions": list(conditions), "rows": rows, "n_pass": int(n_pass), "n_shown": len(rows)}


# --------------------------------------------------------------------------
# Page data
# --------------------------------------------------------------------------

def prepare_data(bigwigs: Dict[str, Dict[str, str]], conditions: Sequence[str], gene_table: pd.DataFrame,
                 gtf: Optional[str] = None, scale_factors: Optional[str] = None, loci_entries: Sequence[str] = (),
                 deg_tables: Optional[Dict[str, str]] = None, top_degs: int = 0, regions_bed: Optional[str] = None,
                 flank: int = 5000, bin_size: int = 100, max_loci: int = 40,
                 metagene_groups: Sequence[str] = ("protein_coding", "lncRNA"), metagene_flank: int = 3000,
                 metagene_min_length: int = 1000, metagene_bed: Optional[str] = None,
                 counts_forward: Optional[str] = None, counts_reverse: Optional[str] = None,
                 counts_matrix: Optional[str] = None, samples: Optional[str] = None, min_count: float = 10,
                 min_samples: int = 2, strand_min_reads: int = 50, threads: int = 1,
                 banner: str = "", title: str = "TTpipe stranded signal") -> Dict[str, Any]:
    """Everything the page shows, as a JSON-serialisable dict (no HTML)."""
    pyBigWig = _import_pybigwig()
    for c in conditions:
        for s in ("forward", "reverse"):
            if s not in bigwigs.get(c, {}):
                raise ValueError(f"no {s} bigWig given for condition {c!r}")
    loci, not_found, notes = select_loci(loci_entries, gene_table, deg_tables, top_degs, regions_bed, max_loci)
    for nf in not_found:
        log(f"NOT FOUND: {nf['entry']}: {nf['reason']}")
    for n in notes:
        log(n)
    opened = {c: {s: pyBigWig.open(bigwigs[c][s]) for s in ("forward", "reverse")} for c in conditions}
    try:
        t0 = time.time()
        records, failed = build_locus_tracks(loci, opened, conditions, flank, bin_size, gtf, gene_table)
        log(f"{len(records)} locus track(s) drawn ({time.time() - t0:.1f} s)")
    finally:
        for c in opened:
            for b in opened[c].values():
                b.close()
    for f in failed:
        log(f"NOT DRAWN: {f['entry']}: {f['reason']}")
    not_found = not_found + failed

    expressed = expressed_gene_ids(counts_matrix, min_count, min_samples)
    if counts_matrix and expressed is not None:
        log(f"expressed genes (count >= {min_count} in >= {min_samples} samples): {len(expressed)}")
    bed_regions = read_bed(metagene_bed) if metagene_bed and os.path.exists(metagene_bed) else None
    mg_notes: List[str] = []
    if metagene_bed and not bed_regions:
        mg_notes.append(f"metagene_bed {metagene_bed!r} is missing or has no usable lines")
    metagene = build_metagene(bigwigs, conditions, gene_table, list(metagene_groups), metagene_flank,
                              metagene_min_length, expressed, bed_regions, threads)
    metagene["notes"] = mg_notes + metagene["notes"]
    for n in metagene["notes"]:
        log(n)

    strand_block = None
    st_note = ""
    if counts_forward and counts_reverse and os.path.exists(counts_forward) and os.path.exists(counts_reverse):
        sc = {}
        if samples and os.path.exists(samples):
            sm = pd.read_csv(samples, sep="\t", dtype=str)
            if {"sample", "condition"} <= set(sm.columns):          # samples_deseq2.txt
                sc = dict(zip(sm["sample"], sm["condition"]))
            elif {"Name", "Condition"} <= set(sm.columns):          # the samplesheet itself
                sc = dict(zip(sm["Name"], sm["Condition"]))
            else:
                raise ValueError(f"{samples}: need columns sample, condition (samples_deseq2.txt) or Name, Condition (samplesheet)")
        stdf, n_pass = strand_table(read_matrix(counts_forward), read_matrix(counts_reverse), gene_table, sc,
                                    conditions, strand_min_reads)
        strand_block = strand_table_block(stdf, conditions, n_pass)
        log(f"strand table: {n_pass} genes with >= {strand_min_reads} reads, {strand_block['n_shown']} embedded")
    else:
        st_note = "counts_forward / counts_reverse were not found: no strand table"

    sfs = read_scale_factors(scale_factors)
    return {
        "meta": {"title": title, "banner": banner or "", "generated": _dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
                 "conditions": list(conditions), "flank": int(flank), "bin_size": int(bin_size),
                 "normalization": NORMALIZATION, "scale_factors": {c: sfs.get(c) for c in conditions},
                 "not_found": not_found, "notes": notes, "strand_note": st_note,
                 "strand_min_reads": int(strand_min_reads)},
        "loci": records,
        "metagene": metagene,
        "strand_table": strand_block,
    }


# --------------------------------------------------------------------------
# Optional static figures (matplotlib)
# --------------------------------------------------------------------------

def _import_matplotlib():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise SystemExit("track_static_formats was set but matplotlib is not installed in this environment "
                         f"(conda: matplotlib-base; envs/ttpipe_plots.yml has it): {e}")
    return plt


def safe_name(label: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", label).strip("_") or "locus"


def write_static_figures(data: Dict[str, Any], outdir: str, formats: Sequence[str]) -> List[str]:
    """One matplotlib figure per locus (same layout as the page: + strand above, - strand mirrored below, one
    panel per condition on a shared x axis and, here, a shared y scale, then gene models)."""
    bad = [f for f in formats if f not in ALLOWED_STATIC]
    if bad:
        raise SystemExit(f"track_static_formats: unknown format(s) {bad}; allowed: {list(ALLOWED_STATIC)}")
    plt = _import_matplotlib()
    from matplotlib.patches import Rectangle
    os.makedirs(outdir, exist_ok=True)
    conds = data["meta"]["conditions"]
    plus_col, minus_col = "#0072B2", "#D55E00"
    written, used = [], set()
    for loc in data["loci"]:
        nrow = max([m["row"] for m in loc["models"]] + [0]) + 1
        ratios = [3] * len(conds) + [0.6 + 0.55 * nrow]
        fig, axes = plt.subplots(len(conds) + 1, 1, sharex=True, figsize=(8, 1.8 * len(conds) + 0.5 * nrow + 1.2),
                                 gridspec_kw={"height_ratios": ratios})
        x0 = loc["win_start"] - 1
        edges = x0 + loc["bin"] * np.arange(loc["n"] + 1)
        fmax = max(max((v or 0) for v in loc["sig"][c]["f"]) for c in conds)
        rmax = max(max((v or 0) for v in loc["sig"][c]["r"]) for c in conds)
        top, bot = max(fmax * 1.08, 1e-9), max(rmax * 1.08, 1e-9)
        for ax, c in zip(axes, conds):
            f = np.array([v or 0 for v in loc["sig"][c]["f"]], float)
            r = np.array([v or 0 for v in loc["sig"][c]["r"]], float)
            ax.stairs(f, edges, fill=True, color=plus_col, linewidth=0)
            ax.stairs(-r, edges, fill=True, color=minus_col, linewidth=0, baseline=0)
            ax.axhline(0, color="#444", linewidth=0.6)
            ax.set_ylim(-bot, top)
            ticks = [t for t in ax.get_yticks() if -bot <= t <= top]
            ax.set_yticks(ticks)
            ax.set_yticklabels([f"{abs(t):g}" for t in ticks], fontsize=7)
            ax.text(0.005, 0.93, c, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")
            ax.text(0.995, 0.93, "+ strand RNA", transform=ax.transAxes, fontsize=7, color=plus_col, ha="right", va="top")
            ax.text(0.995, 0.07, "- strand RNA", transform=ax.transAxes, fontsize=7, color=minus_col, ha="right", va="bottom")
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            ax.set_ylabel("scaled signal", fontsize=7)
        g = axes[-1]
        g.set_ylim(-nrow + 0.4, 0.9)
        g.set_yticks([])
        for sp in ("top", "right", "left"):
            g.spines[sp].set_visible(False)
        for m in loc["models"]:
            y = -m["row"]
            gs, ge = max(m["start"] - 1, edges[0]), min(m["end"], edges[-1])
            g.plot([gs, ge], [y, y], color="#333", linewidth=0.8, zorder=1)
            for s, e in m["exons"]:
                g.add_patch(Rectangle((s - 1, y - 0.22), e - s + 1, 0.44, color="#333", zorder=2, linewidth=0))
            arrow = ">" if m["strand"] == "+" else ("<" if m["strand"] == "-" else "")
            g.text(gs, y + 0.3, f"{m['name']} {arrow}", fontsize=7, va="bottom")
        g.set_xlim(edges[0], edges[-1])
        g.set_xlabel(f"{loc['chrom']} (bp)", fontsize=8)
        g.ticklabel_format(axis="x", style="plain", useOffset=False)
        g.tick_params(axis="x", labelsize=7)
        fig.suptitle(f"{loc['label']}  {loc['chrom']}:{loc['start']}-{loc['end']} ({loc['strand']})", fontsize=10)
        fig.text(0.01, 0.005, "y axis: scaled coverage (library-size scale factors already applied in the bigWigs)", fontsize=5.5, color="#666")
        fig.tight_layout(rect=(0, 0.015, 1, 0.97))
        base = safe_name(loc["label"])
        k = 2
        while base in used:
            base = f"{safe_name(loc['label'])}_{k}"
            k += 1
        used.add(base)
        for fmt in formats:
            p = os.path.join(outdir, f"{base}.{fmt}")
            fig.savefig(p, dpi=200 if fmt == "png" else None)
            written.append(p)
        plt.close(fig)
    return written


# --------------------------------------------------------------------------
# HTML
# --------------------------------------------------------------------------

EXTRA_CSS = """
.pill{display:inline-block;background:var(--card);border:1px solid var(--border);border-radius:10px;padding:0 8px;font-size:.8rem;color:var(--muted);margin-right:4px}
#locusplot,#mgplot{min-height:420px}
.legendline{display:flex;gap:16px;flex-wrap:wrap;font-size:.88rem;margin:4px 0}
.sw{display:inline-block;width:22px;height:0;border-top:3px solid;vertical-align:middle;margin-right:4px}
.nfbox{background:var(--warnbg);color:var(--warn);border:1px solid var(--warn);border-radius:6px;padding:8px 12px;margin:12px 0}
.nfbox ul{margin:4px 0 0 18px;padding:0}
"""

JS = r"""
(function(){
'use strict';
var D = JSON.parse(document.getElementById('ttdata').textContent);
var $ = function(id){ return document.getElementById(id); };
var PLUS = '#0072B2', MINUS = '#D55E00', CPAL = ['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9','#7a5195','#555555'];
var C = D.meta.conditions;
var S = {li:0, freey:false, visible:[], mg:0, mview:'tss', mstat:'mean', tf:{text:'', type:'', min:0, anti:0}, tsort:{key:'total', asc:false}, tshown:200, trows:[]};

function esc(s){ return String(s).replace(/[&<>"']/g, function(c){ return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]; }); }
function fmt(v){ var a = Math.abs(v); if(a === 0) return '0'; if(a >= 1000 || a < 0.01) return v.toExponential(1).replace('e+', 'e'); return String(+v.toPrecision(3)); }
function fire(el, type){ el.dispatchEvent(new Event(type, {bubbles:true})); }
function gridCol(){ return 'rgba(128,128,128,0.25)'; }
function fgc(){ try { return getComputedStyle(document.body).color; } catch(e){ return '#222'; } }
function mx(a){ var m = 0; for(var i = 0; i < a.length; i++){ if(a[i] !== null && a[i] > m) m = a[i]; } return m; }
function niceStep(max){ var raw = max / 3, p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p; return (f < 1.5 ? 1 : (f < 3.5 ? 2 : (f < 7.5 ? 5 : 10))) * p; }
/* y axis with the strand mirrored: tick labels show absolute values, + above 0 and - below 0 */
function yAxisFor(fmax, rmax){
  var top = Math.max(fmax * 1.08, rmax * 0.15, 1e-9), bot = Math.max(rmax * 1.08, fmax * 0.15, 1e-9);
  var step = niceStep(Math.max(fmax, rmax, 1e-9)), vals = [0], text = ['0'], v;
  for(v = step; v <= fmax * 1.001; v += step){ vals.push(v); text.push(fmt(v)); }
  for(v = step; v <= rmax * 1.001; v += step){ vals.push(-v); text.push(fmt(v)); }
  if(rmax > 0 && vals.every(function(t){ return t >= 0; })){ vals.push(-rmax); text.push(fmt(rmax)); }
  if(fmax > 0 && vals.every(function(t){ return t <= 0; })){ vals.push(fmax); text.push(fmt(fmax)); }
  return {range:[-bot, top], tickvals:vals, ticktext:text};
}

/* ---------------- locus tracks ---------------- */
function locusTitle(L){ return L.label + '  ' + L.chrom + ':' + L.start.toLocaleString('en-US') + '-' + L.end.toLocaleString('en-US') + ' (' + L.strand + ')'; }

function drawLocus(){
  var L = D.loci[S.li];
  if(!L){ $('locusplot').innerHTML = ''; return Promise.resolve(); }
  var nC = C.length, rows = 0, i, k;
  L.models.forEach(function(m){ rows = Math.max(rows, m.row + 1); });
  rows = Math.max(rows, 1);
  var sigH = 140, geneH = 34 * rows + 26, top = 56, bot = 62, gap = 10;
  var H = top + bot + nC * sigH + (nC - 1) * gap + gap + geneH;
  function dom(y0, y1){ return [1 - (y1 - top) / (H - top - bot), 1 - (y0 - top) / (H - top - bot)]; }
  var fmaxAll = 0, rmaxAll = 0;
  C.forEach(function(c){ fmaxAll = Math.max(fmaxAll, mx(L.sig[c].f)); rmaxAll = Math.max(rmaxAll, mx(L.sig[c].r)); });
  var shared = yAxisFor(fmaxAll, rmaxAll);
  var layout = {height:H, margin:{l:78, r:20, t:top, b:bot}, showlegend:false, barmode:'overlay', bargap:0, paper_bgcolor:'rgba(0,0,0,0)', plot_bgcolor:'rgba(0,0,0,0)',
                font:{color:fgc(), size:12}, shapes:[], annotations:[], hovermode:'closest'};
  var data = [], xs = [], cur = top;
  for(i = 0; i < L.n; i++) xs.push(L.win_start - 1 + (i + 0.5) * L.bin);
  var xr = [L.win_start - 1, L.win_start - 1 + L.n * L.bin];
  var geneAx = nC + 1;                     // x/y axis number of the gene panel
  var visible = [];
  C.forEach(function(c, ci){
    var n = ci + 1, d = dom(cur, cur + sigH), fa = L.sig[c];
    var ya = S.freey ? yAxisFor(mx(fa.f), mx(fa.r)) : shared;
    layout['yaxis' + (n === 1 ? '' : n)] = {domain:d, range:ya.range, tickvals:ya.tickvals, ticktext:ya.ticktext, zeroline:true, zerolinecolor:fgc(), zerolinewidth:1,
      showgrid:true, gridcolor:gridCol(), anchor:'x' + n, tickfont:{size:10}, fixedrange:false};
    layout['xaxis' + (n === 1 ? '' : n)] = {domain:[0, 1], range:xr, anchor:'y' + n, matches:'x' + geneAx, showticklabels:false, showgrid:false, zeroline:false};
    var neg = fa.r.map(function(v){ return v === null ? null : -v; });
    data.push({type:'bar', x:xs, y:fa.f, width:L.bin, marker:{color:PLUS, line:{width:0}}, name:c + ' + strand RNA', xaxis:'x' + n, yaxis:'y' + n,
               hovertemplate:'%{x:,.0f}<br>+ strand RNA: %{y:.3g}<extra>' + esc(c) + '</extra>'});
    data.push({type:'bar', x:xs, y:neg, width:L.bin, marker:{color:MINUS, line:{width:0}}, name:c + ' - strand RNA', xaxis:'x' + n, yaxis:'y' + n,
               customdata:fa.r, hovertemplate:'%{x:,.0f}<br>- strand RNA: %{customdata:.3g}<extra>' + esc(c) + '</extra>'});
    var ydom = d;
    layout.annotations.push({xref:'paper', yref:'paper', x:0.005, y:ydom[1], xanchor:'left', yanchor:'top', text:'<b>' + esc(c) + '</b>', showarrow:false, font:{size:13}});
    layout.annotations.push({xref:'paper', yref:'paper', x:0.995, y:ydom[1], xanchor:'right', yanchor:'top', text:'+ strand RNA', showarrow:false, font:{size:11, color:PLUS}});
    layout.annotations.push({xref:'paper', yref:'paper', x:0.995, y:ydom[0], xanchor:'right', yanchor:'bottom', text:'- strand RNA', showarrow:false, font:{size:11, color:MINUS}});
    visible.push(c);
    cur += sigH + gap;
  });
  /* gene panel */
  var gd = dom(cur, cur + geneH);
  var gn = geneAx;
  layout['yaxis' + gn] = {domain:gd, range:[-(rows - 1) - 0.65, 0.95], showgrid:false, zeroline:false, showticklabels:false, anchor:'x' + gn, fixedrange:true};
  layout['xaxis' + (gn === 1 ? '' : gn)] = {domain:[0, 1], range:xr, anchor:'y' + gn, title:{text:L.chrom + ' (bp)', standoff:6}, tickformat:',d', showgrid:false, zeroline:false, ticks:'outside'};
  var ar = {x:[], y:[], sym:[], col:[]}, hov = {x:[], y:[], t:[]};
  L.models.forEach(function(m){
    var y = -m.row, gs = Math.max(m.start - 1, xr[0]), ge = Math.min(m.end, xr[1]);
    layout.shapes.push({type:'line', xref:'x' + gn, yref:'y' + gn, x0:gs, x1:ge, y0:y - 0.1, y1:y - 0.1, line:{color:fgc(), width:1}});
    m.exons.forEach(function(e){
      layout.shapes.push({type:'rect', xref:'x' + gn, yref:'y' + gn, x0:e[0] - 1, x1:e[1], y0:y - 0.3, y1:y + 0.08, fillcolor:fgc(), line:{width:0}});
    });
    layout.annotations.push({xref:'x' + gn, yref:'y' + gn, x:gs, y:y + 0.14, xanchor:'left', yanchor:'bottom', text:esc(m.name) + (m.strand === '+' ? ' &#9656;' : (m.strand === '-' ? ' &#9666;' : '')), showarrow:false, font:{size:11}});
    var step = (xr[1] - xr[0]) / 28;
    for(var p = gs + step / 2; p < ge; p += step){
      var inEx = m.exons.some(function(e){ return p >= e[0] - 1 && p <= e[1]; });
      if(m.strand === '+' || m.strand === '-'){ ar.x.push(p); ar.y.push(y - 0.1); ar.sym.push(m.strand === '+' ? 'triangle-right' : 'triangle-left'); ar.col.push(inEx ? '#ffffff' : fgc()); }
    }
    hov.x.push((gs + ge) / 2); hov.y.push(y - 0.1); hov.t.push(esc(m.name) + ' (' + esc(m.id) + ') ' + m.strand + ' ' + L.chrom + ':' + m.start + '-' + m.end);
  });
  if(ar.x.length) data.push({type:'scatter', mode:'markers', x:ar.x, y:ar.y, xaxis:'x' + gn, yaxis:'y' + gn, marker:{symbol:ar.sym, size:6, color:ar.col, line:{width:0}}, hoverinfo:'skip', showlegend:false});
  if(hov.x.length) data.push({type:'scatter', mode:'markers', x:hov.x, y:hov.y, xaxis:'x' + gn, yaxis:'y' + gn, marker:{size:14, opacity:0}, text:hov.t, hovertemplate:'%{text}<extra></extra>', showlegend:false});
  if(L.gene_id){
    layout.shapes.push({type:'rect', xref:'x' + gn, yref:'paper', x0:L.start - 1, x1:L.end, y0:0, y1:1, fillcolor:'rgba(128,128,128,0.12)', line:{width:0}, layer:'below'});
  }
  layout.annotations.push({xref:'paper', yref:'paper', x:-0.062, y:0.5, text:'scaled coverage (+ above, - below the axis)' + (S.freey ? ' - free y per condition' : ' - same y scale in all panels'),
                           textangle:-90, showarrow:false, font:{size:11}, xanchor:'center'});
  layout.annotations.push({xref:'paper', yref:'paper', x:0, y:1 + 34 / (H - top - bot), xanchor:'left', yanchor:'bottom', showarrow:false, font:{size:14},
                           text:'<b>' + esc(locusTitle(L)) + '</b>'});
  S.layout = layout;
  var cfg = {responsive:true, displaylogo:false};
  $('locinfo').innerHTML = '<span class="pill">' + esc(L.chrom) + ':' + L.win_start.toLocaleString('en-US') + '-' + L.win_end.toLocaleString('en-US') + '</span>' +
     '<span class="pill">' + L.bin + ' bp bins' + (L.bin_raised ? ' (raised from ' + D.meta.bin_size + ' bp: wide window)' : '') + '</span>' +
     L.sources.map(function(s){ return '<span class="pill">' + esc(s) + '</span>'; }).join('') +
     (L.models_hidden ? '<span class="pill">' + L.models_hidden + ' more gene(s) not drawn</span>' : '') +
     (L.models.length ? '' : '<span class="pill">no annotated gene in this window</span>');
  return Plotly.react('locusplot', data, layout, cfg);
}

function fillLoci(){
  var sel = $('locus'), q = S.q || '', keep = [];
  sel.innerHTML = '';
  D.loci.forEach(function(L, i){
    var hay = (L.label + ' ' + L.gene_id + ' ' + L.chrom + ':' + L.start + '-' + L.end + ' ' + L.sources.join(' ')).toLowerCase();
    if(!q || hay.indexOf(q) >= 0){ keep.push(i); var o = document.createElement('option'); o.value = i; o.textContent = L.label + '  (' + L.chrom + ':' + L.start + '-' + L.end + ')'; sel.appendChild(o); }
  });
  S.visible = keep;
  $('lsearchmsg').textContent = keep.length + ' of ' + D.loci.length + ' loci';
  if(keep.length){ if(keep.indexOf(S.li) < 0) S.li = keep[0]; sel.value = String(S.li); }
}

/* ---------------- metagene ---------------- */
function drawMeta(){
  var G = D.metagene.groups[S.mg];
  if(!G){ return Promise.resolve(); }
  var V = G[S.mview];
  if(!V){ $('mgnote').textContent = G.body_note || 'This view is not available for this group.'; $('mgnote').hidden = false; $('mgplot').hidden = true; return Promise.resolve(); }
  $('mgnote').hidden = true; $('mgplot').hidden = false;
  var data = [], fm = 0, rm = 0, ci = 0, st = S.mstat, sense = G.oriented ? 'sense' : '+ strand RNA', anti = G.oriented ? 'antisense' : '- strand RNA';
  C.forEach(function(c, i){
    var p = V.profiles[c], s = p['sense_' + st], a = p['anti_' + st], col = CPAL[i % CPAL.length];
    fm = Math.max(fm, mx(s)); rm = Math.max(rm, mx(a));
    data.push({type:'scatter', mode:'lines', x:V.x, y:s, name:c + ' ' + sense, line:{color:col, width:2}, legendgroup:c, hovertemplate:'%{y:.3g}<extra>' + esc(c) + ' ' + sense + '</extra>'});
    data.push({type:'scatter', mode:'lines', x:V.x, y:a.map(function(v){ return v === null ? null : -v; }), name:c + ' ' + anti, line:{color:col, width:2, dash:'dash'},
               legendgroup:c, customdata:a, hovertemplate:'%{customdata:.3g}<extra>' + esc(c) + ' ' + anti + '</extra>'});
  });
  var ya = yAxisFor(fm, rm);
  var title = G.label + ' - ' + (S.mview === 'tss' ? (G.kind === 'bed' ? 'centred on the region centre' : 'centred on the TSS') : 'gene body scaled to ' + D.metagene.body_bins + ' bins, with flanks') + ' - ' + st;
  var layout = {height:480, margin:{l:78, r:20, t:56, b:60}, paper_bgcolor:'rgba(0,0,0,0)', plot_bgcolor:'rgba(0,0,0,0)', font:{color:fgc(), size:12},
    xaxis:{range:(S.mview === 'tss' ? [V.ticks.vals[0], V.ticks.vals[V.ticks.vals.length - 1]] : [-0.5, V.x.length - 0.5]), tickvals:V.ticks.vals, ticktext:V.ticks.text, showgrid:false, zeroline:false, title:{text:S.mview === 'tss' ? (G.kind === 'bed' ? 'position relative to the region centre (5\' to 3\' of the region if it has a strand)' : 'position relative to the TSS (5\' on the left)') : '5\' to 3\' (gene orientation)'}},
    yaxis:{range:ya.range, tickvals:ya.tickvals, ticktext:ya.ticktext, zeroline:true, zerolinecolor:fgc(), gridcolor:gridCol(), title:{text:'scaled coverage per bin'}},
    legend:{orientation:'h', y:-0.22}, title:{text:esc(title), font:{size:14}, x:0.01, xanchor:'left'}, hovermode:'x unified', shapes:[],
    annotations:[{xref:'paper', yref:'y', x:0.995, xanchor:'right', y:ya.range[1], yanchor:'top', text:esc(sense) + ' (solid, above)', showarrow:false, font:{size:11}},
                 {xref:'paper', yref:'y', x:0.995, xanchor:'right', y:ya.range[0], yanchor:'bottom', text:esc(anti) + ' (dashed, mirrored below)', showarrow:false, font:{size:11}}]};
  if(S.mview === 'body'){ layout.shapes.push({type:'rect', xref:'x', yref:'paper', x0:V.ticks.vals[1], x1:V.ticks.vals[2], y0:0, y1:1, fillcolor:'rgba(128,128,128,0.10)', line:{width:0}, layer:'below'}); }
  else { layout.shapes.push({type:'line', xref:'x', yref:'paper', x0:0, x1:0, y0:0, y1:1, line:{color:'rgba(128,128,128,0.6)', dash:'dot', width:1}}); }
  var n = S.mview === 'tss' ? G.n_tss : G.n_body;
  $('mginfo').textContent = n + ' ' + (G.kind === 'bed' ? 'regions' : 'genes') + ' in this view' + (G.kind === 'gene' ? ' (expressed genes of type ' + G.name + ': ' + G.n_expressed + ')' : '') + (G.kind === 'bed' && !G.oriented ? '; no strand column: strands shown as + / -' : '');
  return Plotly.react('mgplot', data, layout, {responsive:true, displaylogo:false});
}

function fillMetaControls(){
  var sel = $('mgroup'); sel.innerHTML = '';
  D.metagene.groups.forEach(function(g, i){ var o = document.createElement('option'); o.value = i; o.textContent = g.label; sel.appendChild(o); });
  var any = D.metagene.groups.length > 0;
  $('mgcontrols').hidden = !any; $('mgnone').hidden = any;
  $('mgnotes').innerHTML = D.metagene.notes.map(function(n){ return '<li>' + esc(n) + '</li>'; }).join('');
  $('mgnotes').parentNode.parentNode.hidden = D.metagene.notes.length === 0;
}

/* ---------------- strand table ---------------- */
var TCOLS = [['name','Gene','l'],['type','Type','l'],['loc','Locus','l'],['strand','Strand',''],['total','Total reads','']];
function tval(r, key){
  var T = D.strand_table;
  if(key === 'name') return r[1].toLowerCase();
  if(key === 'type') return r[2].toLowerCase();
  if(key === 'loc') return r[3] + ':' + r[4];
  if(key === 'strand') return r[6];
  if(key === 'total') return r[7];
  if(key === 'pooled'){ var s = 0, a = 0; r[8].forEach(function(x){ s += x[0]; a += x[1]; }); return s + a > 0 ? a / (s + a) : null; }
  var m = key.match(/^c(\d+)$/); if(m) return r[8][+m[1]][2];
  return null;
}
function pooled(r){ return tval(r, 'pooled'); }
function tableRows(){
  var T = D.strand_table, f = S.tf, q = f.text.toLowerCase(), out = [];
  T.rows.forEach(function(r, i){
    if(q && (r[0] + ' ' + r[1] + ' ' + r[2]).toLowerCase().indexOf(q) < 0) return;
    if(f.type && r[2] !== f.type) return;
    if(r[7] < f.min) return;
    if(f.anti > 0){ var ok = r[8].some(function(x){ return x[2] !== null && x[2] >= f.anti; }); if(!ok) return; }
    out.push(i);
  });
  var key = S.tsort.key, d = S.tsort.asc ? 1 : -1;
  out.sort(function(a, b){
    var x = tval(T.rows[a], key), y = tval(T.rows[b], key);
    if(x === null && y === null) return a - b; if(x === null) return 1; if(y === null) return -1;
    if(x < y) return -d; if(x > y) return d; return a - b;
  });
  return out;
}
function fmtFrac(v){ return v === null ? 'NA' : (v < 0.01 && v > 0 ? v.toFixed(4) : v.toFixed(3)); }
function drawTable(){
  var T = D.strand_table; if(!T) return;
  var rows = tableRows(); S.trows = rows;
  var h = '<thead><tr>';
  function th(key, label, cls){ var a = S.tsort.key === key ? (S.tsort.asc ? 'ascending' : 'descending') : 'none'; return '<th class="' + cls + '" aria-sort="' + a + '"><button type="button" data-sort="' + key + '">' + esc(label) + '</button></th>'; }
  TCOLS.forEach(function(c){ h += th(c[0], c[1], c[2]); });
  h += th('pooled', 'Antisense fraction (all samples)', '');
  T.conditions.forEach(function(c, i){ h += th('c' + i, 'Antisense fraction ' + c, ''); });
  h += '</tr></thead><tbody>';
  rows.slice(0, S.tshown).forEach(function(i){
    var r = T.rows[i];
    h += '<tr data-gi="' + i + '"><td class="l">' + esc(r[1]) + '</td><td class="l">' + esc(r[2]) + '</td><td class="l">' + esc(r[3]) + ':' + r[4] + '-' + r[5] + '</td><td>' + esc(r[6]) + '</td><td>' + r[7] + '</td><td>' + fmtFrac(pooled(r)) + '</td>';
    r[8].forEach(function(x){ h += '<td title="sense ' + x[0] + ', antisense ' + x[1] + '">' + fmtFrac(x[2]) + '</td>'; });
    h += '</tr>';
  });
  $('stbl').innerHTML = h + '</tbody>';
  $('stinfo').textContent = rows.length + ' of ' + T.n_shown + ' embedded genes match (' + T.n_pass + ' genes have at least ' + D.meta.strand_min_reads + ' reads in total; the ' + T.n_shown + ' with the most reads are embedded)';
  $('stmore').hidden = rows.length <= S.tshown;
  Array.prototype.forEach.call(document.querySelectorAll('#stbl th button'), function(b){
    b.onclick = function(){ var k = b.getAttribute('data-sort'); if(S.tsort.key === k) S.tsort.asc = !S.tsort.asc; else { S.tsort.key = k; S.tsort.asc = (k === 'name' || k === 'type' || k === 'loc' || k === 'strand'); } S.tshown = 200; drawTable(); };
  });
}
function csvEsc(v){ var s = String(v); return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; }
function currentCsv(){
  var T = D.strand_table, head = ['gene_id','gene_name','gene_type','chr','start','end','strand','total_reads','antisense_fraction_all_samples'];
  T.conditions.forEach(function(c){ head.push(c + '_sense_reads', c + '_antisense_reads', c + '_antisense_fraction'); });
  var lines = [head.map(csvEsc).join(',')];
  S.trows.forEach(function(i){
    var r = T.rows[i], p = pooled(r), row = [r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], p === null ? '' : p.toFixed(4)];
    r[8].forEach(function(x){ row.push(x[0], x[1], x[2] === null ? '' : x[2]); });
    lines.push(row.map(csvEsc).join(','));
  });
  return lines.join('\n') + '\n';
}

/* ---------------- wiring ---------------- */
function showTab(id){
  Array.prototype.forEach.call(document.querySelectorAll('.tabs button'), function(b){ b.setAttribute('aria-selected', b.getAttribute('data-panel') === id ? 'true' : 'false'); });
  Array.prototype.forEach.call(document.querySelectorAll('.panel'), function(p){ p.hidden = p.id !== id; });
  var w = window.Plotly;
  if(id === 'p-loci' && w) { drawLocus(); }
  if(id === 'p-meta' && w) { drawMeta(); }
}
function csvDownload(){
  var blob = new Blob([currentCsv()], {type:'text/csv'}), a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = 'TTpipe_antisense_fraction.csv'; document.body.appendChild(a); a.click(); document.body.removeChild(a);
}

function init(){
  $('nfbox').hidden = D.meta.not_found.length === 0;
  $('nflist').innerHTML = D.meta.not_found.map(function(n){ return '<li><code>' + esc(n.entry) + '</code>: ' + esc(n.reason) + '</li>'; }).join('');
  $('locnotes').innerHTML = D.meta.notes.map(function(n){ return '<li>' + esc(n) + '</li>'; }).join('');
  $('locnotes').parentNode.parentNode.hidden = D.meta.notes.length === 0;
  Array.prototype.forEach.call(document.querySelectorAll('.tabs button'), function(b){ b.onclick = function(){ showTab(b.getAttribute('data-panel')); }; });
  fillLoci(); fillMetaControls();
  $('locus').onchange = function(){ S.li = +$('locus').value; drawLocus(); };
  $('lprev').onclick = function(){ var k = S.visible.indexOf(S.li); if(k > 0){ S.li = S.visible[k - 1]; $('locus').value = String(S.li); drawLocus(); } };
  $('lnext').onclick = function(){ var k = S.visible.indexOf(S.li); if(k >= 0 && k < S.visible.length - 1){ S.li = S.visible[k + 1]; $('locus').value = String(S.li); drawLocus(); } };
  $('freey').onchange = function(){ S.freey = $('freey').checked; drawLocus(); };
  $('lsearch').oninput = function(){ S.q = $('lsearch').value.trim().toLowerCase(); fillLoci(); drawLocus(); };
  $('mgroup').onchange = function(){ S.mg = +$('mgroup').value; drawMeta(); };
  Array.prototype.forEach.call(document.querySelectorAll('input[name=mview]'), function(r){ r.onchange = function(){ if(r.checked){ S.mview = r.value; drawMeta(); } }; });
  Array.prototype.forEach.call(document.querySelectorAll('input[name=mstat]'), function(r){ r.onchange = function(){ if(r.checked){ S.mstat = r.value; drawMeta(); } }; });
  var T = D.strand_table;
  $('p-table-none').hidden = !!T;
  $('p-table-body').hidden = !T;
  if(T){
    var types = {}; T.rows.forEach(function(r){ types[r[2]] = (types[r[2]] || 0) + 1; });
    var ks = Object.keys(types).sort(function(a, b){ return types[b] - types[a]; });
    $('ttype').innerHTML = '<option value="">all gene types</option>' + ks.map(function(k){ return '<option value="' + esc(k) + '">' + esc(k || 'unknown') + ' (' + types[k] + ')</option>'; }).join('');
    $('ttext').oninput = function(){ S.tf.text = $('ttext').value; S.tshown = 200; drawTable(); };
    $('ttype').onchange = function(){ S.tf.type = $('ttype').value; S.tshown = 200; drawTable(); };
    $('tmin').onchange = function(){ S.tf.min = +$('tmin').value || 0; S.tshown = 200; drawTable(); };
    $('tanti').onchange = function(){ S.tf.anti = +$('tanti').value || 0; S.tshown = 200; drawTable(); };
    $('stmore').onclick = function(){ S.tshown += 200; drawTable(); };
    $('stcsv').onclick = csvDownload;
    drawTable();
  }
  if(!D.loci.length){ $('locusplot').hidden = true; $('noloci').hidden = false; }
  return drawLocus();
}

function showErrors(){
  var errs = window.__ttErrors || [], box = $('errbox');
  if(errs.length){ box.hidden = false; box.textContent = 'JavaScript error(s): ' + errs.join(' | '); }
}

/* ---------------- optional self-test: open the file with #selftest ---------------- */
function selfTest(){
  var out = [], fails = 0, chain = Promise.resolve();
  function ok(name, cond, extra){ if(!cond) fails++; out.push((cond ? 'PASS ' : 'FAIL ') + name + (extra ? ' ' + extra : '')); }
  function step(f){ chain = chain.then(f); }
  function wait(ms){ return new Promise(function(res){ setTimeout(res, ms); }); }
  function same(a, b){ return a.length === b.length && a.every(function(v, i){ return Math.abs(v - b[i]) < 1e-9; }); }
  function barTraces(){ return $('locusplot').data.filter(function(t){ return t.type === 'bar'; }); }
  step(function(){
    var P = $('locusplot'), L = D.loci[0], bars = barTraces();
    ok('locus options', $('locus').options.length === D.loci.length, String($('locus').options.length));
    ok('two bar traces per condition', bars.length === 2 * C.length, String(bars.length));
    ok('+ strand traces are all >= 0 and come from the forward data', bars.filter(function(t, i){ return i % 2 === 0; }).every(function(t, ci){ return t.y.every(function(v){ return v >= 0; }) && same(t.y, L.sig[C[ci]].f); }));
    ok('- strand traces are mirrored below zero and come from the reverse data', bars.filter(function(t, i){ return i % 2 === 1; }).every(function(t, ci){ return t.y.every(function(v){ return v <= 0; }) && same(t.y.map(function(v){ return -v; }), L.sig[C[ci]].r); }));
    var labs = P.layout.annotations.map(function(a){ return a.text; });
    ok('strand labels present', labs.indexOf('+ strand RNA') >= 0 && labs.indexOf('- strand RNA') >= 0);
    var r1 = P.layout.yaxis.range, rl = P.layout['yaxis' + C.length].range;
    ok('y scale shared across conditions', same(r1, rl), JSON.stringify(r1));
    ok('x axes are linked (shared x)', P.layout['xaxis2'] ? P.layout.xaxis2.matches === 'x' + (C.length + 1) : true);
    ok('gene models drawn', L.models.length === 0 || P.layout.shapes.some(function(s){ return s.type === 'rect' && s.yref === 'y' + (C.length + 1); }), L.models.length + ' model(s)');
  });
  step(function(){
    var fy = $('freey'); fy.checked = true; fire(fy, 'change');
    return wait(120).then(function(){
      var P = $('locusplot'), a = P.layout.yaxis.range, b = P.layout['yaxis' + C.length].range;
      ok('free y toggle gives per-condition ranges', C.length < 2 || !same(a, b), JSON.stringify(a) + ' vs ' + JSON.stringify(b));
      fy.checked = false; fire(fy, 'change'); return wait(120);
    }).then(function(){
      var P = $('locusplot');
      ok('shared y restored', same(P.layout.yaxis.range, P.layout['yaxis' + C.length].range));
    });
  });
  step(function(){
    if(D.loci.length < 2){ ok('second locus (only one locus in this page)', true); return; }
    var sel = $('locus'); sel.value = '1'; fire(sel, 'change');
    return wait(150).then(function(){
      var P = $('locusplot'), L = D.loci[1];
      ok('locus dropdown redraws', P.layout.xaxis.range[0] === L.win_start - 1 && P.layout.annotations.some(function(a){ return a.text.indexOf(L.label) >= 0; }), L.label);
      ok('bins match the second locus', barTraces()[0].x.length === L.n, String(L.n));
    });
  });
  step(function(){
    if(D.loci.length < 2){ return; }
    var q = D.loci[D.loci.length - 1].label, s = $('lsearch');
    s.value = q.toLowerCase(); fire(s, 'input');
    var n = $('locus').options.length;
    ok('search keeps matching loci only', n >= 1 && n < D.loci.length || D.loci.length === 1, n + ' of ' + D.loci.length);
    return wait(120).then(function(){
      ok('search selects a matching locus', D.loci[S.li].label.toLowerCase().indexOf(q.toLowerCase()) >= 0);
      s.value = 'zzzz_no_such_locus'; fire(s, 'input');
      ok('search with no match empties the list', $('locus').options.length === 0, $('lsearchmsg').textContent);
      s.value = ''; fire(s, 'input'); return wait(100);
    }).then(function(){ ok('search cleared', $('locus').options.length === D.loci.length); });
  });
  step(function(){
    var b = document.querySelector('.tabs button[data-panel=p-meta]'); b.click();
    return wait(200).then(function(){
      var G = D.metagene.groups;
      if(!G.length){ ok('metagene note shown when no group qualifies', !$('mgnone').hidden); return; }
      var P = $('mgplot');
      ok('metagene drawn: sense and antisense line per condition', P.data.length === 2 * C.length, String(P.data.length));
      ok('antisense is mirrored below zero and dashed', P.data[1].y.every(function(v){ return v <= 0; }) && P.data[1].line.dash === 'dash' && P.data[0].line.dash !== 'dash');
      var sm = G[0].tss.profiles[C[0]].sense_mean;
      ok('sense mean line equals the embedded data', same(P.data[0].y, sm));
      var radio = document.querySelector('input[name=mstat][value=median]'); radio.checked = true; fire(radio, 'change');
      return wait(120).then(function(){
        ok('median toggle switches the line', same($('mgplot').data[0].y, G[0].tss.profiles[C[0]].sense_median));
        var m = document.querySelector('input[name=mstat][value=mean]'); m.checked = true; fire(m, 'change'); return wait(100);
      }).then(function(){
        ok('mean toggle switches back', same($('mgplot').data[0].y, sm));
        if(G[0].body){
          var v = document.querySelector('input[name=mview][value=body]'); v.checked = true; fire(v, 'change');
          return wait(120).then(function(){
            ok('gene-body view drawn', same($('mgplot').data[0].y, G[0].body.profiles[C[0]].sense_mean) && $('mgplot').data[0].x.length === G[0].body.x.length);
            var t = document.querySelector('input[name=mview][value=tss]'); t.checked = true; fire(t, 'change'); return wait(100);
          });
        }
      }).then(function(){
        if(G.length > 1){
          var sel = $('mgroup'); sel.value = '1'; fire(sel, 'change');
          return wait(150).then(function(){
            ok('metagene group switch', $('mgplot').layout.title.text.indexOf(G[1].label) >= 0 && same($('mgplot').data[0].y, G[1].tss.profiles[C[0]].sense_mean), G[1].label);
          });
        } else { ok('single metagene group in this page', true); }
      });
    });
  });
  step(function(){
    var T = D.strand_table;
    if(!T){ ok('no strand table in this page', !$('p-table-none').hidden); return; }
    document.querySelector('.tabs button[data-panel=p-table]').click();
    var n0 = document.querySelectorAll('#stbl tbody tr').length;
    ok('strand table rows drawn', n0 > 0, String(n0));
    var b = document.querySelector('#stbl th button[data-sort=pooled]'); b.click();
    var rows = S.trows, good = true, d = S.tsort.asc ? 1 : -1;
    for(var i = 1; i < rows.length; i++){ var x = tval(T.rows[rows[i - 1]], 'pooled'), y = tval(T.rows[rows[i]], 'pooled'); if(x !== null && y !== null && (y - x) * d < 0) good = false; }
    ok('table sorts by antisense fraction (' + (S.tsort.asc ? 'ascending' : 'descending') + ')', S.tsort.key === 'pooled' && good, rows.length + ' rows');
    document.querySelector('#stbl th button[data-sort=pooled]').click();
    var asc2 = S.tsort.asc; rows = S.trows; good = true; d = asc2 ? 1 : -1;
    for(var j = 1; j < rows.length; j++){ var x2 = tval(T.rows[rows[j - 1]], 'pooled'), y2 = tval(T.rows[rows[j]], 'pooled'); if(x2 !== null && y2 !== null && (y2 - x2) * d < 0) good = false; }
    ok('second click reverses the sort', good);
    var tt = $('ttext'); tt.value = T.rows[0][1]; fire(tt, 'input');
    ok('table text filter', S.trows.length >= 1 && S.trows.length < T.rows.length || T.rows.length === 1, S.trows.length + ' rows for ' + T.rows[0][1]);
    var csv = currentCsv().trim().split('\n');
    ok('csv rows match the shown rows', csv.length - 1 === S.trows.length && csv[0].indexOf('gene_id,gene_name,gene_type') === 0, (csv.length - 1) + ' rows');
    tt.value = 'zzzz_none'; fire(tt, 'input');
    ok('filter with no match', S.trows.length === 0);
    tt.value = ''; fire(tt, 'input');
    var mn = $('tmin'); mn.value = '1e12'; fire(mn, 'change');
    ok('min total reads filter', S.trows.length === 0);
    mn.value = '0'; fire(mn, 'change');
    ok('csv text is available for all rows again', currentCsv().trim().split('\n').length - 1 === S.trows.length && S.trows.length === T.rows.length);
  });
  step(function(){ ok('no script errors', (window.__ttErrors || []).length === 0, (window.__ttErrors || []).join('|')); });
  return chain.catch(function(e){ fails++; out.push('FAIL exception ' + e + (e && e.stack ? ' ' + e.stack.split('\n')[1] : '')); }).then(function(){
    out.push('SELFTEST DONE failures=' + fails);
    var pre = $('selftest'); pre.hidden = false; pre.textContent = out.join('\n');
  });
}

try {
  init().then(function(){
    document.body.setAttribute('data-ready', '1');
    showErrors();
    if(location.hash === '#selftest') return selfTest();
  }).catch(function(e){ (window.__ttErrors = window.__ttErrors || []).push(String(e) + (e && e.stack ? ' ' + e.stack.split('\n')[1] : '')); showErrors(); });
} catch(e){ (window.__ttErrors = window.__ttErrors || []).push(String(e)); showErrors(); }
window.addEventListener('error', showErrors);
})();
"""

BODY = """
<div class="wrap">
<header><h1>__TITLE__</h1><div class="sub">__SUB__</div></header>
__BANNER__
<div id="errbox" class="err" role="alert" hidden></div>
<div class="note"><b>Normalization.</b> __NORM__ __SF__</div>
<div id="nfbox" class="nfbox" role="note" hidden><b>Not found / not drawn</b><ul id="nflist"></ul></div>
<details open><summary>How to read this page</summary>
<div class="note">
<p><b>Locus tracks.</b> One panel per condition. The RNA made from the + genomic strand (<span style="color:#0072B2"><b>+ strand RNA</b></span>, the
<code>forward</code> bigWig) is drawn above the axis and the RNA made from the - strand (<span style="color:#D55E00"><b>- strand RNA</b></span>, the
<code>reverse</code> bigWig) is drawn mirrored below it. Both bigWigs hold positive numbers; the minus sign only places the - strand below the axis, and tick
labels and hover text show the real value. A gene on the + strand has its sense signal above the axis, a gene on the - strand below it. Gene models underneath
are the union of the annotated exons of each gene (introns as a line, arrow = strand); the grey band is the gene itself. All panels share the x axis (zoom in
one, all follow) and, by default, the same y scale, so conditions can be compared; tick "free y" to scale each panel on its own.</p>
<p><b>Metagene.</b> Mean (or median) signal per bin over many genes, relative to gene orientation (5' on the left). Solid line = sense, dashed line below the
axis = antisense (the opposite strand), one colour per condition. <b>Strand table.</b> Antisense fraction = antisense reads / (sense + antisense reads) from the
HTSeq counts. TTpipe does not call enhancer RNAs or any other transcript class: a BED file you provide only tells it where to look.</p>
</div></details>
<div class="tabs" role="tablist">
<button type="button" role="tab" aria-selected="true" data-panel="p-loci">Locus tracks</button>
<button type="button" role="tab" aria-selected="false" data-panel="p-meta">Metagene</button>
<button type="button" role="tab" aria-selected="false" data-panel="p-table">Strand table</button>
</div>
<section id="p-loci" class="panel">
<details hidden><summary>Notes</summary><div class="note"><ul id="locnotes"></ul></div></details>
<div id="noloci" class="note" hidden>No locus was selected or could be drawn. Set <code>track_loci</code>, <code>track_regions_bed</code> or <code>track_top_degs</code> in the config and rerun this rule.</div>
<div class="controls">
<div class="ctl"><label for="locus">Locus</label><div class="row"><button type="button" class="btn" id="lprev" aria-label="previous locus">&lt;</button>
<select id="locus"></select><button type="button" class="btn" id="lnext" aria-label="next locus">&gt;</button></div></div>
<div class="ctl"><label for="lsearch">Search the selected loci (name, gene_id, coordinates, source)</label><div class="row"><input id="lsearch" type="text" autocomplete="off" placeholder="e.g. Gapdh"><span id="lsearchmsg" class="sub"></span></div></div>
<div class="ctl"><span>Y scale</span><div class="row"><label><input id="freey" type="checkbox"> free y (each condition on its own scale)</label></div></div>
</div>
<div id="locinfo" class="sub"></div>
<div id="locusplot"></div>
</section>
<section id="p-meta" class="panel" hidden>
<details hidden><summary>Notes</summary><div class="note"><ul id="mgnotes"></ul></div></details>
<div id="mgnone" class="note" hidden>No metagene group qualified (see the notes above): groups need at least 20 expressed genes, and the bigWigs must cover the genes' windows.</div>
<div id="mgcontrols" class="controls">
<div class="ctl"><label for="mgroup">Gene group (GENCODE-style gene_type, expressed genes)</label><select id="mgroup"></select></div>
<div class="ctl"><span>View</span><div class="row"><label><input type="radio" name="mview" value="tss" checked> TSS-centred</label>
<label><input type="radio" name="mview" value="body"> gene body scaled</label></div></div>
<div class="ctl"><span>Statistic across genes</span><div class="row"><label><input type="radio" name="mstat" value="mean" checked> mean</label>
<label><input type="radio" name="mstat" value="median"> median</label></div></div>
</div>
<div id="mginfo" class="sub"></div>
<div id="mgnote" class="note" hidden></div>
<div id="mgplot"></div>
</section>
<section id="p-table" class="panel" hidden>
<div id="p-table-none" class="note" hidden>__STNOTE__</div>
<div id="p-table-body">
<p class="note">Antisense fraction = antisense reads / (sense + antisense reads), counting reads on the whole gene (introns included) from the strand-split HTSeq counts, summed over the samples of a condition.
For a + gene sense = forward counts, for a - gene sense = reverse counts. Only genes with at least __MINREADS__ reads over all samples are listed; the top __TOPN__ by total reads are embedded.
A value near 0 means almost all reads are on the gene's own strand; antisense transcription, overlapping genes on the other strand, or a strandedness mismatch raise it.
Hover a cell for the read counts. NA = no reads in that condition.</p>
<div class="controls">
<div class="ctl"><label for="ttext">Filter (name, gene_id, type)</label><input id="ttext" type="text"></div>
<div class="ctl"><label for="ttype">Gene type</label><select id="ttype"></select></div>
<div class="ctl"><label for="tmin">Min total reads</label><input id="tmin" type="number" min="0" step="any" value="0"></div>
<div class="ctl"><label for="tanti">Min antisense fraction (in any condition)</label><input id="tanti" type="number" min="0" max="1" step="0.05" value="0"></div>
<div class="ctl"><span>&nbsp;</span><button type="button" class="btn" id="stcsv">Download CSV of the shown rows</button></div>
</div>
<div class="sub" id="stinfo"></div>
<div class="tablewrap"><table id="stbl"></table></div>
<p><button type="button" class="btn" id="stmore" hidden>Show 200 more rows</button> <span class="sub">The CSV contains every row matching the filters, not only the rows displayed.</span></p>
</div>
</section>
<pre id="selftest" hidden></pre>
</div>
"""


def _esc(x: Any) -> str:
    import html
    return html.escape(str(x), quote=True)


def build_html(data: Dict[str, Any], plotly_js: Optional[str] = None) -> str:
    """The complete page. `plotly_js` defaults to the plotly.js bundled with the installed plotly package."""
    if plotly_js is None:
        from plotly.offline import get_plotlyjs
        plotly_js = get_plotlyjs()
    m = data["meta"]
    n_mg = len(data["metagene"]["groups"])
    sub = (f"{len(m['conditions'])} condition(s) ({_esc(', '.join(m['conditions']))}) &middot; {len(data['loci'])} locus track(s) "
           f"&middot; {n_mg} metagene group(s) &middot; window = locus &plusmn; {m['flank']} bp, {m['bin_size']} bp bins "
           f"&middot; generated {_esc(m['generated'])}")
    banner = f'<div class="banner" role="note">{_esc(m["banner"])}</div>' if m.get("banner") else ""
    sf = ", ".join(f"{_esc(c)} = {_esc(round(v, 4))}" for c, v in m["scale_factors"].items() if v is not None)
    sf_txt = f"Library-size scale factors already applied: {sf}." if sf else ""
    st = data.get("strand_table")
    body = (BODY.replace("__TITLE__", _esc(m["title"])).replace("__SUB__", sub).replace("__BANNER__", banner)
            .replace("__NORM__", _esc(m["normalization"])).replace("__SF__", sf_txt)
            .replace("__STNOTE__", _esc(m.get("strand_note") or "No strand table."))
            .replace("__MINREADS__", str(m["strand_min_reads"])).replace("__TOPN__", "5000"))
    head_js = "window.__ttErrors=[];window.addEventListener('error',function(e){window.__ttErrors.push((e.message||'error')+(e.lineno?' @'+e.lineno:''));});"
    css = BASE_CSS + EXTRA_CSS
    return ("<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">"
            f"<meta name=\"color-scheme\" content=\"light dark\"><title>{_esc(m['title'])}</title>\n<style>{css}</style>\n"
            f"<script>{head_js}</script>\n<script>{plotly_js}</script></head>\n<body>{body}\n"
            f"<script id=\"ttdata\" type=\"application/json\">{to_json(data)}</script>\n<script>{JS}</script>\n</body></html>\n")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def parse_bigwigs(specs: Sequence[str]) -> Dict[str, Dict[str, str]]:
    out: Dict[str, Dict[str, str]] = {}
    for s in specs or []:
        parts = s.split("=", 2)
        if len(parts) != 3 or parts[1] not in ("forward", "reverse") or not parts[0] or not parts[2]:
            raise ValueError(f"--bigwig expects CONDITION=forward|reverse=PATH, got {s!r}")
        out.setdefault(parts[0], {})[parts[1]] = parts[2]
    return out


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Stranded signal plots (locus tracks, metagene, strand table).")
    ap.add_argument("--bigwig", nargs="+", action="extend", default=None, help="CONDITION=forward|reverse=PATH, one or more (may be repeated)")
    ap.add_argument("--condition-order", default="", help="comma-separated; default: order of --bigwig")
    ap.add_argument("--scale-factors", default=None)
    ap.add_argument("--gene-info", required=True)
    ap.add_argument("--gtf", default=None)
    ap.add_argument("--counts-forward", default=None)
    ap.add_argument("--counts-reverse", default=None)
    ap.add_argument("--counts-matrix", default=None, help="count matrix used for the expressed-gene filter")
    ap.add_argument("--samples", default=None, help="samples_deseq2.txt (sample, condition) or the samplesheet (Name, Condition)")
    ap.add_argument("--loci", nargs="*", default=[], help="gene symbols, gene_ids or chr:start-end")
    ap.add_argument("--degs", nargs="*", default=[], help="DESeq2 result csv files (degs_<comparison>.csv)")
    ap.add_argument("--top-degs", type=int, default=0)
    ap.add_argument("--regions-bed", default="")
    ap.add_argument("--flank", type=int, default=5000)
    ap.add_argument("--bin-size", type=int, default=100)
    ap.add_argument("--max-loci", type=int, default=40)
    ap.add_argument("--static-formats", nargs="*", default=[], help="any of svg pdf png (checked in the code)")
    ap.add_argument("--static-dir", default="")
    ap.add_argument("--metagene-groups", nargs="*", default=["protein_coding", "lncRNA"])
    ap.add_argument("--metagene-flank", type=int, default=3000)
    ap.add_argument("--metagene-min-length", type=int, default=1000)
    ap.add_argument("--metagene-bed", default="")
    ap.add_argument("--strand-table-min-reads", type=int, default=50)
    ap.add_argument("--min-count", type=float, default=10)
    ap.add_argument("--min-samples", type=int, default=2)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--banner", default="")
    ap.add_argument("--title", default="TTpipe stranded signal")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    bad = [f for f in a.static_formats if f not in ALLOWED_STATIC]   # checked here, not with argparse's choices option: that fails with nargs="*" defaults on Python 3.11
    if bad:
        ap.error(f"unknown --static-formats {bad}; choose from {list(ALLOWED_STATIC)}")
    if a.static_formats:
        _import_matplotlib()                      # fail early with a clear message
    if a.bin_size < 1 or a.flank < 0 or a.metagene_flank < TSS_BIN_BP:
        ap.error("--bin-size must be >= 1, --flank >= 0 and --metagene-flank >= %d" % TSS_BIN_BP)
    bws = parse_bigwigs(a.bigwig)
    if not bws:
        ap.error("no --bigwig given")
    conds = [c for c in a.condition_order.split(",") if c] or list(bws)
    conds = [c for c in conds if c in bws] + [c for c in bws if c not in conds]
    gt = read_gene_table(a.gene_info)
    degs = {contrast_name_from_path(p): p for p in a.degs}
    data = prepare_data(bws, conds, gt, gtf=a.gtf or None, scale_factors=a.scale_factors, loci_entries=a.loci,
                        deg_tables=degs, top_degs=a.top_degs, regions_bed=a.regions_bed or None, flank=a.flank,
                        bin_size=a.bin_size, max_loci=a.max_loci, metagene_groups=a.metagene_groups,
                        metagene_flank=a.metagene_flank, metagene_min_length=a.metagene_min_length,
                        metagene_bed=a.metagene_bed or None, counts_forward=a.counts_forward,
                        counts_reverse=a.counts_reverse, counts_matrix=a.counts_matrix, samples=a.samples,
                        min_count=a.min_count, min_samples=a.min_samples, strand_min_reads=a.strand_table_min_reads,
                        threads=a.threads, banner=a.banner, title=a.title)
    doc = build_html(data)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(doc)
    mb = os.path.getsize(a.out) / 1e6
    log(f"wrote {a.out} ({mb:.1f} MB; {len(data['loci'])} locus track(s), {len(data['metagene']['groups'])} metagene "
        f"group(s), strand table: {'yes' if data['strand_table'] else 'no'}, {len(data['meta']['not_found'])} not found)")
    if mb > SIZE_WARN_MB:
        log(f"WARNING: the file is {mb:.1f} MB (> {SIZE_WARN_MB:.0f} MB); fewer loci (track_max_loci) or larger bins make it smaller.")
    if a.static_formats:
        outdir = a.static_dir or os.path.join(os.path.dirname(os.path.abspath(a.out)), "stranded_signal")
        files = write_static_figures(data, outdir, a.static_formats)
        log(f"wrote {len(files)} static figure file(s) to {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
