#!/usr/bin/env python3
"""Static, self-contained TTpipe QC report (one HTML file, plotly.js inlined; works offline).

Usage:
    python scripts/qc_report.py --outdir OUTDIR --samplesheet samplesheet.txt \
        --config config/TTpipeConfig.yaml [more.yaml ...] --out OUTDIR/qc/TTpipe_qc_report.html \
        [--sample-qc OUTDIR/qc/sample_qc.tsv] [--counts MATRIX] [--run-processing 1] [--run-analysis 1] \
        [--banner "text shown at the top"]

Files read (all optional; paths follow workflows/rules/common.smk under OUTDIR):
    qc/sample_qc.tsv                 SAMPLE_QC, made by scripts/qc_summary.py (the per-sample metrics)
    counts/counts_deseq2.txt or --counts   gene x sample matrix: correlation heatmap and (fallback) PCA
    analysis/deseq2/vst.tsv          VST_COUNTS: PCA when the analysis ran
    analysis/deseq2/deg_summary.tsv  DEG_SUMMARY: comparison, n_up, n_down, n_tested
    tracks/scale_factors.tsv         SCALE_FACTORS: condition, mapped_pairs, scale_factor, source
    tracks/**/*.bw                   bigWig files (only listed, no genome plots)
Thresholds come from `qc_thresholds` in the config (same semantics as MicroCPipe: higher_is_better fails below
`fail` and warns below `warn`; lower_is_better the reverse; a value exactly on an edge is not penalised).

The report never fabricates numbers: a section without data is replaced by a note, and every skipped or
partial section is listed with its reason near the top. Things the report does not show (spike-in fraction,
which the optional spike-in module writes to spike/spike_counts.tsv; labeled-vs-total enrichment;
synthesis/degradation rates) are listed in a "Not assessed" box, not drawn.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import glob
import html
import math
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qc_summary as qs  # noqa: E402  (shared parsers: same directory)

DEFAULT_THRESHOLDS: Dict[str, Dict[str, Any]] = {
    "input_pairs_million": {"direction": "higher_is_better", "warn": 20, "fail": 5},
    "unique_pct": {"direction": "higher_is_better", "warn": 70, "fail": 50},
    "multi_pct": {"direction": "lower_is_better", "warn": 15, "fail": 30},
    "too_short_pct": {"direction": "lower_is_better", "warn": 20, "fail": 40},
    "duplicate_pct": {"direction": "lower_is_better", "warn": 40, "fail": 50},
    "sense_pct": {"direction": "higher_is_better", "warn": 90, "fail": 75},
    "rrna_pct": {"direction": "lower_is_better", "warn": 5, "fail": 15},
    "replicate_correlation": {"direction": "higher_is_better", "warn": 0.95, "fail": 0.85},
}

METRIC_INFO: Dict[str, Dict[str, str]] = {
    "input_pairs_million": {"label": "Input pairs", "unit": "M", "digits": "1",
                            "desc": "STAR 'Number of input reads' (read pairs), in millions."},
    "unique_pct": {"label": "Unique", "unit": "%", "digits": "1",
                   "desc": "STAR uniquely mapped reads, % of input pairs."},
    "multi_pct": {"label": "Multi", "unit": "%", "digits": "1",
                  "desc": "STAR reads mapped to multiple + too many loci, % of input pairs."},
    "too_short_pct": {"label": "Too short", "unit": "%", "digits": "1",
                      "desc": "STAR reads unmapped because the alignment was too short, % of input pairs."},
    "duplicate_pct": {"label": "Duplicates", "unit": "%", "digits": "1",
                      "desc": "Picard PERCENT_DUPLICATION x 100, as reported by Picard MarkDuplicates."},
    "sense_pct": {"label": "Sense", "unit": "%", "digits": "1",
                  "desc": "% of counted reads on the annotated gene strand (HTSeq strand-split counts; STAR "
                          "ReadsPerGene when those are absent). Checks the 'strandedness' setting."},
    "rrna_pct": {"label": "rRNA", "unit": "%", "digits": "2",
                 "desc": "% of counted reads in genes annotated rRNA / Mt_rRNA. mm10 references annotate few "
                         "rRNA loci, so this underestimates real rRNA contamination."},
    "replicate_correlation": {"label": "Replicate r", "unit": "", "digits": "3",
                              "desc": "Mean Pearson r of log2(CPM+1) gene counts against the other replicates of "
                                      "the same condition."},
}
EXTRA_INFO = {
    "other_unmapped_pct": ("Other unmapped", "%", "STAR unmapped: too many mismatches + other, % of input pairs (not thresholded)."),
    "final_pct_of_input": ("Final / input", "%", "Read pairs in the final deduplicated BAM / STAR input pairs (not thresholded)."),
}
STATUS_RANK = {"pass": 0, "warn": 1, "fail": 2, "na": -1}
STATUS_SYMBOL = {"pass": "✓", "warn": "!", "fail": "✕", "na": "–"}
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#999999"]

NOT_ASSESSED = [
    ("Spike-in fraction", "not shown in this report. With the optional spike-in module on (spike_in.enabled), the spike-in read pairs and percentage per sample are in spike/spike_counts.tsv and per spike-in contig in spike/spike_contig_counts.tsv."),
    ("Labeled-vs-total enrichment", "this version has no total / input RNA samples to compare the labeled (4sU) RNA against."),
    ("Synthesis and degradation rates", "rates need total and labeled RNA plus spike-in normalisation, none of which are modelled here."),
]

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------


def esc(x: Any) -> str:
    return html.escape(str(x), quote=True)


def isnum(v: Any) -> bool:
    return v is not None and not (isinstance(v, float) and math.isnan(v)) and not pd.isna(v)


def fmt(v: Any, digits: int = 1) -> str:
    return f"{float(v):,.{digits}f}" if isnum(v) else "n/a"


def classify(value: Optional[float], rule: Dict[str, Any]) -> str:
    """'pass' / 'warn' / 'fail' / 'na'. higher_is_better: < fail fails, < warn warns; lower_is_better: > fail, > warn."""
    if not isnum(value):
        return "na"
    warn, fail = rule.get("warn"), rule.get("fail")
    if rule.get("direction", "higher_is_better") == "higher_is_better":
        if fail is not None and value < fail:
            return "fail"
        if warn is not None and value < warn:
            return "warn"
    else:
        if fail is not None and value > fail:
            return "fail"
        if warn is not None and value > warn:
            return "warn"
    return "pass"


def worst_status(statuses) -> str:
    best = "na"
    for s in statuses:
        if STATUS_RANK[s] > STATUS_RANK[best]:
            best = s
    return best


def classify_sample(metrics: Dict[str, Optional[float]], thresholds: Dict[str, Dict[str, Any]]) -> Dict[str, str]:
    out = {m: classify(metrics.get(m), r) for m, r in thresholds.items()}
    out["overall"] = worst_status(out.values())
    return out


def load_thresholds(config: dict) -> Dict[str, Dict[str, Any]]:
    user = (config or {}).get("qc_thresholds") or {}
    thr = {k: dict(v) for k, v in user.items() if isinstance(v, dict)}
    return thr or {k: dict(v) for k, v in DEFAULT_THRESHOLDS.items()}


def _read_tsv(path: str) -> Optional[pd.DataFrame]:
    if not os.path.exists(path):
        return None
    try:
        return pd.read_csv(path, sep="\t")
    except Exception as exc:
        print(f"qc_report: could not read {path}: {exc}", file=sys.stderr)
        return None


_AXIS = dict(gridcolor="rgba(128,128,128,0.25)", linecolor="rgba(128,128,128,0.6)",
             zerolinecolor="rgba(128,128,128,0.4)", automargin=True)


def _layout(fig: go.Figure, height: int = 420, **kw: Any) -> go.Figure:
    fig.update_layout(height=height, margin=dict(l=60, r=20, t=30, b=50), paper_bgcolor="rgba(0,0,0,0)",
                      plot_bgcolor="rgba(0,0,0,0)", font=dict(family="system-ui, sans-serif", size=13),
                      legend=dict(orientation="h", yanchor="top", y=-0.25, x=0), hoverlabel=dict(namelength=-1), **kw)
    fig.update_xaxes(**_AXIS)
    fig.update_yaxes(**_AXIS)
    return fig


def _fig_html(fig: go.Figure, key: str) -> str:
    return '<div class="plot">' + fig.to_html(full_html=False, include_plotlyjs=False, div_id=f"plot-{key}",
                                              config={"responsive": True, "displaylogo": False}) + "</div>"


def cond_colours(conditions: List[str]) -> Dict[str, str]:
    return {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(conditions)}


# --------------------------------------------------------------------------
# PCA
# --------------------------------------------------------------------------


def compute_pca(mat: pd.DataFrame, n_top: int = 1000) -> Optional[Tuple[pd.DataFrame, np.ndarray]]:
    """PCA of samples (columns of mat; already on a log-like scale). Uses the n_top most variable genes,
    centred per gene. Returns (coords DataFrame indexed by sample with PC1.., explained-variance fractions) or
    None if fewer than 3 samples or no variable genes."""
    if mat.shape[1] < 3:
        return None
    x = mat.to_numpy(dtype=float).T  # samples x genes
    var = x.var(axis=0)
    keep = np.where(var > 0)[0]
    if keep.size < 2:
        return None
    keep = keep[np.argsort(var[keep])[::-1][:n_top]]
    x = x[:, keep] - x[:, keep].mean(axis=0)
    u, s, _vt = np.linalg.svd(x, full_matrices=False)
    npc = min(len(s), 4)
    coords = pd.DataFrame(u[:, :npc] * s[:npc], index=mat.columns, columns=[f"PC{i + 1}" for i in range(npc)])
    expl = (s ** 2) / (s ** 2).sum()
    return coords, expl[:npc]


# --------------------------------------------------------------------------
# Sections. Each returns HTML and records skips in `skips` as (section, reason)
# --------------------------------------------------------------------------

Skips = List[Tuple[str, str]]


def _note(text: str) -> str:
    return f'<p class="note">{esc(text)}</p>'


def _badge(status: str, text: Optional[str] = None, cls: str = "") -> str:
    word = text or {"pass": "PASS", "warn": "WARN", "fail": "FAIL", "na": "N/A"}[status]
    return (f'<span class="badge {status} {cls}"><span class="sym" aria-hidden="true">{STATUS_SYMBOL[status]}</span>'
            f"{esc(word)}</span>")


def _threshold_text(rule: Dict[str, Any]) -> str:
    higher = rule.get("direction", "higher_is_better") == "higher_is_better"
    return f"pass {'>=' if higher else '<='} {rule.get('warn')}, fail {'<' if higher else '>'} {rule.get('fail')}"


def _why_missing(what: str, ctx: Dict[str, Any]) -> str:
    """Reason text for absent per-sample processing statistics."""
    if ctx["run_processing"] is False:
        return f"{what} not available: run_processing is false (analysis-only run)."
    return f"{what} not found under {ctx['outdir']}."


def has(df: pd.DataFrame, col: str) -> bool:
    return col in df.columns and df[col].notna().any()


def render_header(ctx: Dict[str, Any], qc: pd.DataFrame) -> str:
    cfg = ctx["config"]
    conds = ctx["conditions"]
    per = {c: int((ctx["sheet"]["Condition"] == c).sum()) for c in conds}
    rows = [
        ("Samples", f"{len(ctx['sheet'])} in {len(conds)} condition(s): "
                    + ", ".join(f"{c} (n={n})" for c, n in per.items())),
        ("Mode", f"processing: {'yes' if ctx['run_processing'] else 'no'}; analysis: "
                 f"{'yes' if ctx['run_analysis'] else 'no'}"),
        ("Genome", cfg.get("genome", "n/a")),
        ("Library strandedness", cfg.get("strandedness", "dUTP")),
        ("Count strand for DESeq2", cfg.get("count_strand", "sum")),
        ("Trimming (fastp)", "yes" if cfg.get("trim") else "no"),
        ("MAPQ filter / dedup", f"MAPQ >= {cfg.get('min_mapq', 'n/a')}; duplicates "
                                f"{'removed' if cfg.get('remove_duplicates', True) else 'marked only'}"),
        ("DESeq2", f"reference {cfg.get('reference_condition', 'n/a')}; design {cfg.get('design', 'n/a')}; "
                   f"alpha {cfg.get('alpha', 'n/a')}"),
        ("Output directory", ctx["outdir"]),
    ]
    dl = "".join(f"<dt>{esc(k)}</dt><dd>{esc(v)}</dd>" for k, v in rows)
    return f'<section id="run"><h2>Run</h2><dl>{dl}</dl></section>'


def render_summary(qc: pd.DataFrame, thresholds, statuses, skips: Skips, ctx) -> str:
    counts = {k: sum(1 for st in statuses.values() if st["overall"] == k) for k in ("pass", "warn", "fail", "na")}
    cards = "".join(f'<div class="card"><b>{counts[k]}</b>{_badge(k)}</div>' for k in ("pass", "warn", "fail", "na")
                    if counts[k] or k != "na")
    cards = f'<div class="cards">{cards}<div class="card"><b>{len(qc)}</b>sample(s)</div></div>'
    btn = '<th scope="col"{t}><button type="button">{l}</button></th>'
    cols = list(thresholds)
    heads = [btn.format(t="", l="Sample"), btn.format(t="", l="Condition"), btn.format(t="", l="Status")]
    for m in cols:
        info = METRIC_INFO.get(m, {"label": m, "unit": "", "desc": "(custom metric)"})
        unit = f" ({info['unit']})" if info["unit"] else ""
        heads.append(btn.format(t=f' title="{esc(info["desc"])} Thresholds: {esc(_threshold_text(thresholds[m]))}"',
                                l=esc(info["label"] + unit)))
    extras = [k for k in EXTRA_INFO if k in qc.columns]
    for k in extras:
        lab, unit, desc = EXTRA_INFO[k]
        heads.append(btn.format(t=f' title="{esc(desc)}"', l=esc(f"{lab} ({unit})")))
    rows = []
    for _, r in qc.iterrows():
        st = statuses[r["sample"]]
        ov = st["overall"]
        cells = [f'<td data-sort="{esc(r["sample"])}"><b>{esc(r["sample"])}</b></td>',
                 f'<td class="l" data-sort="{esc(r["condition"])}">{esc(r["condition"])}</td>',
                 f'<td data-sort="{STATUS_RANK[ov]}">{_badge(ov)}</td>']
        for m in cols:
            v = r.get(m)
            sv = st[m]
            dg = int(METRIC_INFO.get(m, {}).get("digits", 1))
            tag = f'<span class="tag">{STATUS_SYMBOL[sv]} {sv.upper()}</span>' if sv in ("warn", "fail") else ""
            cells.append(f'<td class="cell-{sv}" data-sort="{esc(v) if isnum(v) else ""}">{fmt(v, dg)}{tag}</td>')
        for k in extras:
            v = r.get(k)
            cells.append(f'<td data-sort="{esc(v) if isnum(v) else ""}">{fmt(v, 1)}</td>')
        rows.append("<tr>" + "".join(cells) + "</tr>")
    table = ('<div class="tablewrap"><table class="sortable"><thead><tr>' + "".join(heads) + "</tr></thead><tbody>"
             + "".join(rows) + "</tbody></table></div>"
             '<p class="sub legend">Click a header to sort; hover it for the definition and thresholds. n/a = no data '
             f'for that metric (never estimated). {_badge("pass")} {_badge("warn")} {_badge("fail")} {_badge("na")}</p>')
    return f'<section id="summary"><h2>Per-sample QC</h2>{cards}{table}</section>'


def _hline(fig: go.Figure, rule: Optional[Dict[str, Any]], horizontal: bool = False) -> None:
    if not rule:
        return
    for key, dash in (("warn", "dash"), ("fail", "dot")):
        if rule.get(key) is not None:
            kw = dict(line_dash=dash, line_color="rgba(128,128,128,0.9)", annotation_text=key,
                      annotation_position="top left" if not horizontal else "top right")
            (fig.add_vline if horizontal else fig.add_hline)(**{("x" if horizontal else "y"): rule[key]}, **kw)


def _bar(qc: pd.DataFrame, col: str, ytitle: str, colours: Dict[str, str], rule, key: str, yrange=None) -> str:
    d = qc[qc[col].notna()]
    fig = go.Figure(go.Bar(x=d["sample"], y=d[col], marker_color=[colours.get(c, "#999") for c in d["condition"]],
                           hovertemplate="%{x}<br>%{y:.3f}<extra></extra>"))
    _hline(fig, rule)
    _layout(fig, 380, yaxis_title=ytitle, showlegend=False)
    if yrange:
        fig.update_yaxes(range=yrange)
    return _fig_html(fig, key)


def render_alignment(qc, ctx, thresholds, skips: Skips) -> str:
    cols = [("unique_pct", "Uniquely mapped", "#0072B2"), ("multi_pct", "Multi-mapped", "#E69F00"),
            ("too_short_pct", "Unmapped: too short", "#CC79A7"), ("other_unmapped_pct", "Unmapped: other", "#999999")]
    body = ""
    if not has(qc, "unique_pct"):
        reason = _why_missing("STAR Log.final.out files", ctx)
        skips.append(("Alignment", reason))
        body += _note("Skipped. " + reason)
    else:
        fig = go.Figure()
        for c, lab, colr in cols:
            if has(qc, c):
                fig.add_trace(go.Bar(x=qc["sample"], y=qc[c], name=lab, marker_color=colr,
                                     hovertemplate="%{x}<br>" + lab + ": %{y:.2f}%<extra></extra>"))
        _layout(fig, 420, barmode="stack", yaxis_title="% of input read pairs")
        body += _fig_html(fig, "align")
        if has(qc, "mapq_pairs") or has(qc, "final_pairs"):
            fig2 = go.Figure()
            for c, lab, colr in (("input_pairs", "STAR input", "#999999"), ("mapq_pairs", "after MAPQ filter", "#56B4E9"),
                                 ("final_pairs", "final (deduplicated)", "#0072B2")):
                if has(qc, c):
                    fig2.add_trace(go.Bar(x=qc["sample"], y=qc[c] / 1e6, name=lab, marker_color=colr))
            _layout(fig2, 380, barmode="group", yaxis_title="read pairs (millions)")
            body += "<h3>Read-pair retention</h3>" + _fig_html(fig2, "retention")
        else:
            skips.append(("Read-pair retention", _why_missing("MAPQ / final-pair statistics", ctx)))
    return f'<section id="alignment"><h2>Alignment</h2><p class="sub">STAR statistics, % of input read pairs.</p>{body}</section>'


def render_single_metric(sec_id, title, desc, col, ytitle, qc, ctx, thresholds, colours, skips, extra="", yrange=None) -> str:
    if not has(qc, col):
        reason = _why_missing(f"'{col}' data", ctx)
        skips.append((title, reason))
        return f'<section id="{sec_id}"><h2>{esc(title)}</h2><p class="sub">{desc}</p>{_note("Skipped. " + reason)}</section>'
    return (f'<section id="{sec_id}"><h2>{esc(title)}</h2><p class="sub">{desc}</p>'
            + _bar(qc, col, ytitle, colours, thresholds.get(col), sec_id, yrange) + extra + "</section>")


def render_strand(qc, ctx, thresholds, colours, skips: Skips) -> str:
    cfg = ctx["config"]
    desc = (f"% of counted reads on the annotated gene strand; config strandedness = <code>{esc(cfg.get('strandedness', 'dUTP'))}</code>. "
            "A correctly set library gives a high value (dashed/dotted lines: warn/fail thresholds).")
    extra = ""
    if has(qc, "sense_pct"):
        src = sorted(set(qc["sense_source"].dropna()))
        extra = f'<p class="sub">Source: {esc(", ".join(src))} (htseq = gene-body HTSeq counts on the strand-split BAMs; star = STAR ReadsPerGene, exonic reads only).</p>'
        med = float(qc["sense_pct"].median())
        if med < 25:
            extra += _note(f"Median sense is {med:.1f}%: most reads map to the opposite strand. This looks like the "
                           "other library type than the configured strandedness; check the 'strandedness' setting.")
    return render_single_metric("strand", "Strandedness", desc, "sense_pct", "sense reads (%)", qc, ctx, thresholds,
                                colours, skips, extra, yrange=[0, 100])


def render_correlation(qc, ctx, thresholds, colours, skips: Skips) -> str:
    body = ""
    mat = ctx.get("counts")
    if mat is None:
        reason = "no count matrix (counts/counts_deseq2.txt or analysis_counts) found."
        skips.append(("Sample correlation", reason))
        return f'<section id="correlation"><h2>Sample correlation</h2>{_note("Skipped. " + reason)}</section>'
    order = [s for c in ctx["conditions"] for s in ctx["sheet"].loc[ctx["sheet"]["Condition"] == c, "Name"] if s in mat.columns]
    if len(order) < 2:
        reason = "the count matrix has fewer than 2 columns matching the samplesheet Names."
        skips.append(("Sample correlation", reason))
        return f'<section id="correlation"><h2>Sample correlation</h2>{_note("Skipped. " + reason)}</section>'
    lc = qs.log_cpm(mat[order])
    r = lc.corr(method="pearson")
    fig = go.Figure(go.Heatmap(z=r.to_numpy(), x=order, y=order, colorscale="Viridis",
                               zmin=float(max(0.0, np.nanmin(r.to_numpy()))), zmax=1.0, colorbar=dict(title="Pearson r"),
                               hovertemplate="%{y} vs %{x}<br>r = %{z:.4f}<extra></extra>"))
    _layout(fig, max(420, 40 * len(order) + 200))
    fig.update_yaxes(autorange="reversed")
    body += (f'<p class="sub">Pearson r of log2(CPM+1) across all genes ({mat.shape[0]:,} genes), samples ordered by condition. '
             f'Matrix: <code>{esc(ctx["counts_path"])}</code>.</p>' + _fig_html(fig, "corr"))
    if has(qc, "replicate_correlation"):
        body += ("<h3>Replicate correlation</h3><p class=\"sub\">Mean Pearson r vs the other replicates of the same condition "
                 "(pairs with fewer than 3 expressed genes or constant counts are skipped).</p>"
                 + _bar(qc, "replicate_correlation", "mean replicate r", colours, thresholds.get("replicate_correlation"),
                        "repcorr"))
    else:
        reason = "no condition has two or more replicates in the count matrix."
        skips.append(("Replicate correlation", reason))
        body += _note("Replicate correlation skipped. " + reason)
    return f'<section id="correlation"><h2>Sample correlation</h2>{body}</section>'


def render_pca(ctx, colours, skips: Skips) -> str:
    sheet = ctx["sheet"]
    src, mat = None, None
    vst = ctx["vst"]
    if vst is not None and vst.shape[1] >= 3:
        src, mat = f"VST counts from DESeq2 ({ctx['vst_path']}), most variable genes (up to 1000)", vst
    elif ctx.get("counts") is not None and ctx["counts"].shape[1] >= 3:
        mat = qs.log_cpm(ctx["counts"])
        src = (f"log2(CPM+1) computed here from {ctx['counts_path']} (the DESeq2 VST output "
               f"{'is missing' if ctx['run_analysis'] else 'is not available: run_analysis is false'}), most variable genes (up to 1000)")
    if mat is None:
        reason = "needs at least 3 samples in VST counts or in the count matrix; neither is available."
        skips.append(("PCA", reason))
        return f'<section id="pca"><h2>PCA</h2>{_note("Skipped. " + reason)}</section>'
    cols = [s for s in sheet["Name"] if s in mat.columns]
    res = compute_pca(mat[cols])
    if res is None:
        reason = "too few samples or no variable genes for a PCA."
        skips.append(("PCA", reason))
        return f'<section id="pca"><h2>PCA</h2>{_note("Skipped. " + reason)}</section>'
    coords, expl = res
    cond = dict(zip(sheet["Name"], sheet["Condition"]))
    rep = dict(zip(sheet["Name"], sheet["Replicate"]))
    fig = go.Figure()
    for c in ctx["conditions"]:
        ss = [s for s in cols if cond[s] == c]
        if not ss:
            continue
        fig.add_trace(go.Scatter(x=coords.loc[ss, "PC1"], y=coords.loc[ss, "PC2"], mode="markers+text", name=c, text=ss,
                                 textposition="top center", marker=dict(size=12, color=colours[c]),
                                 customdata=[rep[s] for s in ss],
                                 hovertemplate="%{text}<br>replicate %{customdata}<br>PC1 %{x:.2f}, PC2 %{y:.2f}<extra></extra>"))
    _layout(fig, 480, xaxis_title=f"PC1 ({100 * expl[0]:.1f}%)", yaxis_title=f"PC2 ({100 * expl[1]:.1f}%)")
    return f'<section id="pca"><h2>PCA</h2><p class="sub">Input: {esc(src)}. Samples coloured by condition.</p>{_fig_html(fig, "pca")}</section>'


def render_deg(ctx, skips: Skips) -> str:
    d = ctx["deg"]
    if d is None:
        reason = ("run_analysis is false." if ctx["run_analysis"] is False else f"{ctx['deg_path']} not found.")
        skips.append(("Differential expression", reason))
        return f'<section id="deg"><h2>Differential expression</h2>{_note("Skipped. " + reason)}</section>'
    fig = go.Figure()
    fig.add_trace(go.Bar(x=d["comparison"], y=d["n_up"], name="Up", marker_color="#D55E00"))
    fig.add_trace(go.Bar(x=d["comparison"], y=-d["n_down"], name="Down", marker_color="#0072B2",
                         customdata=d["n_down"], hovertemplate="%{x}<br>down: %{customdata}<extra></extra>"))
    _layout(fig, 380, barmode="relative", yaxis_title="DEGs (down negative)")
    tbl = ('<div class="tablewrap"><table><thead><tr><th>Comparison</th><th>Up</th><th>Down</th><th>Tested</th></tr></thead><tbody>'
           + "".join(f"<tr><td>{esc(r.comparison)}</td><td>{int(r.n_up)}</td><td>{int(r.n_down)}</td>"
                     f"<td>{int(r.n_tested) if isnum(getattr(r, 'n_tested', None)) else 'n/a'}</td></tr>"
                     for r in d.itertuples()) + "</tbody></table></div>")
    return (f'<section id="deg"><h2>Differential expression</h2><p class="sub">Genes passing the configured adjusted '
            f'p-value cut-off (alpha {esc(ctx["config"].get("alpha", "n/a"))}) per comparison, from {esc(ctx["deg_path"])}.</p>'
            + _fig_html(fig, "deg") + tbl + "</section>")


def render_tracks(ctx, skips: Skips) -> str:
    sf = ctx["scale"]
    bws = ctx["bigwigs"]
    if sf is None and not bws:
        reason = ("run_processing is false: no tracks were built." if ctx["run_processing"] is False
                  else "no tracks/scale_factors.tsv or bigWig files found.")
        skips.append(("Tracks", reason))
        return f'<section id="tracks"><h2>Tracks</h2>{_note("Skipped. " + reason)}</section>'
    out = '<p class="sub">Listing only: no genome-browser plots are drawn here.</p>'
    if sf is not None:
        head = "".join(f"<th>{esc(c)}</th>" for c in sf.columns)
        body = "".join("<tr>" + "".join(f"<td class=\"l\">{esc(v)}</td>" for v in row) + "</tr>" for row in sf.itertuples(index=False))
        out += f'<h3>Scale factors</h3><div class="tablewrap"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'
    else:
        skips.append(("Scale factors", "tracks/scale_factors.tsv not found."))
    if bws:
        rows = "".join(f"<tr><td class=\"l\"><code>{esc(os.path.relpath(p, ctx['outdir']))}</code></td>"
                       f"<td>{os.path.getsize(p) / 1e6:,.1f}</td></tr>" for p in bws)
        out += (f'<h3>bigWig files ({len(bws)})</h3><div class="tablewrap"><table><thead><tr><th>File</th><th>Size (MB)</th></tr>'
                f"</thead><tbody>{rows}</tbody></table></div>")
    else:
        skips.append(("bigWig list", "no .bw files found under tracks/."))
        out += _note("No bigWig files found under tracks/.")
    return f'<section id="tracks"><h2>Tracks</h2>{out}</section>'


def render_not_assessed() -> str:
    items = "".join(f"<li><b>{esc(a)}</b>: {esc(b)}</li>" for a, b in NOT_ASSESSED)
    return (f'<section id="notassessed"><h2>Not assessed by this pipeline</h2><div class="note"><ul>{items}</ul>'
            "No panels are shown for these because there are no data to show.</div></section>")


def render_methods(thresholds) -> str:
    items = []
    for m, rule in thresholds.items():
        info = METRIC_INFO.get(m, {"desc": "(custom metric)"})
        items.append(f"<dt><code>{esc(m)}</code></dt><dd>{esc(info['desc'])} <span class=\"sub\">Thresholds: "
                     f"{esc(_threshold_text(rule))} ({esc(rule.get('direction', 'higher_is_better').replace('_', ' '))}).</span></dd>")
    for k, (_l, _u, desc) in EXTRA_INFO.items():
        items.append(f"<dt><code>{esc(k)}</code></dt><dd>{esc(desc)}</dd>")
    return f"""<section id="methods"><h2>Methods &amp; metric definitions</h2>
<p>Status rule: for <em>higher is better</em> metrics a value below <code>fail</code> fails and below <code>warn</code> warns;
for <em>lower is better</em> the reverse; a value exactly on an edge is not penalised. A sample's status is its worst metric.
Thresholds are <code>qc_thresholds</code> from the config; the <code>warn</code> values follow the TT-seq guide's rules of thumb,
the <code>fail</code> values are TTpipe defaults, not literature values. Metrics without data are n/a and do not affect the status.</p>
<dl>{''.join(items)}</dl></section>"""


# --------------------------------------------------------------------------
# Page
# --------------------------------------------------------------------------

CSS = """
:root{--bg:#fff;--fg:#1b1f24;--muted:#5b6570;--card:#f6f8fa;--border:#d0d7de;--accent:#0072B2;
--pass:#007a5a;--warn:#8a5a00;--fail:#b3470a;--na:#5b6570;--passbg:#d9f2ea;--warnbg:#fbeccb;--failbg:#f9dccb;--nabg:#e8e8e8}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#0f1318;--fg:#e6e9ed;--muted:#9aa5b1;--card:#171c23;
--border:#2f3741;--accent:#56B4E9;--pass:#3ed6a5;--warn:#f4c04d;--fail:#ff9a5c;--na:#9aa5b1;
--passbg:#0b3a2d;--warnbg:#4a3a0a;--failbg:#4d230c;--nabg:#2a2f36}}
:root[data-theme="dark"]{--bg:#0f1318;--fg:#e6e9ed;--muted:#9aa5b1;--card:#171c23;--border:#2f3741;--accent:#56B4E9;
--pass:#3ed6a5;--warn:#f4c04d;--fail:#ff9a5c;--na:#9aa5b1;--passbg:#0b3a2d;--warnbg:#4a3a0a;--failbg:#4d230c;--nabg:#2a2f36}
*{box-sizing:border-box}html{scroll-behavior:smooth}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}
.wrap{max-width:1250px;margin:0 auto;padding:0 16px 64px}
header{padding:24px 0 8px;display:flex;flex-wrap:wrap;gap:16px;align-items:center;justify-content:space-between}
h1{font-size:1.6rem;margin:0}h2{font-size:1.25rem;margin:40px 0 8px;border-bottom:1px solid var(--border);padding-bottom:6px}
h3{font-size:1.02rem;margin:22px 0 4px}
.sub{color:var(--muted);font-size:.9rem}
nav{position:sticky;top:0;z-index:5;background:var(--bg);border-bottom:1px solid var(--border);padding:8px 0;display:flex;flex-wrap:wrap;gap:4px 16px;font-size:.9rem}
nav a{color:var(--accent);text-decoration:none}nav a:hover{text-decoration:underline}
.badge{display:inline-flex;align-items:center;gap:6px;border-radius:999px;padding:2px 12px;font-weight:600;font-size:.85rem;border:1px solid currentColor;white-space:nowrap}
.badge .sym{font-weight:800}
.pass{color:var(--pass);background:var(--passbg)}.warn{color:var(--warn);background:var(--warnbg)}
.fail{color:var(--fail);background:var(--failbg)}.na{color:var(--na);background:var(--nabg)}
.overall{font-size:1.1rem;padding:6px 18px}
.cards{display:flex;flex-wrap:wrap;gap:10px;margin:8px 0}
.card{background:var(--card);border:1px solid var(--border);border-radius:8px;padding:8px 14px;min-width:110px}
.card b{display:block;font-size:1.4rem}
.tablewrap{overflow-x:auto;border:1px solid var(--border);border-radius:8px}
table{border-collapse:collapse;width:100%;font-size:.9rem}
th,td{padding:7px 10px;text-align:right;border-bottom:1px solid var(--border);white-space:nowrap}
th:first-child,td:first-child{text-align:left;position:sticky;left:0;background:var(--bg)}
thead th{background:var(--card)}
td.l{text-align:left}
th button{all:unset;cursor:pointer;font-weight:600;display:inline-flex;gap:4px;align-items:center}
th button:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
th[aria-sort="ascending"] button::after{content:"\\25B2";font-size:.7em}
th[aria-sort="descending"] button::after{content:"\\25BC";font-size:.7em}
td.cell-pass{background:var(--passbg)}td.cell-warn{background:var(--warnbg)}td.cell-fail{background:var(--failbg)}
td .tag{font-size:.75em;font-weight:700;margin-left:4px}
tr:last-child td{border-bottom:none}
.note{background:var(--card);border-left:4px solid var(--na);padding:8px 12px;border-radius:4px;color:var(--muted)}
.banner{background:var(--warnbg);color:var(--warn);border:1px solid var(--warn);border-radius:6px;padding:8px 12px;margin:12px 0;font-weight:600}
.legend span{margin-right:10px}
code{background:var(--card);padding:1px 5px;border-radius:4px}
.plot{margin:6px 0}
dl{display:grid;grid-template-columns:max-content 1fr;gap:4px 16px}dt{font-weight:600}dd{margin:0}
@media (max-width:640px){dl{grid-template-columns:1fr}h1{font-size:1.3rem}}
.js-plotly-plot .plotly text{fill:var(--fg)!important}
.js-plotly-plot .plotly .xtick text,.js-plotly-plot .plotly .ytick text,.js-plotly-plot .plotly .annotation-text{fill:var(--muted)!important}
"""

JS = """
(function(){
  document.querySelectorAll('table.sortable').forEach(function(tb){
    var ths = tb.querySelectorAll('thead th');
    ths.forEach(function(th, idx){
      var btn = th.querySelector('button'); if(!btn) return;
      btn.addEventListener('click', function(){
        var asc = th.getAttribute('aria-sort') !== 'ascending';
        ths.forEach(function(o){o.removeAttribute('aria-sort');});
        th.setAttribute('aria-sort', asc ? 'ascending' : 'descending');
        var body = tb.tBodies[0], rows = Array.prototype.slice.call(body.rows);
        rows.sort(function(a, b){
          var x = a.cells[idx].getAttribute('data-sort'), y = b.cells[idx].getAttribute('data-sort');
          var nx = parseFloat(x), ny = parseFloat(y);
          if(isNaN(nx) && !isNaN(ny)) return 1;
          if(!isNaN(nx) && isNaN(ny)) return -1;
          var r = (!isNaN(nx) && !isNaN(ny)) ? nx - ny : String(x).localeCompare(String(y));
          return asc ? r : -r;
        });
        rows.forEach(function(r){body.appendChild(r);});
      });
    });
  });
})();
"""


def load_context(outdir: str, sheet: pd.DataFrame, config: dict, counts_path: Optional[str],
                 run_processing: Optional[bool], run_analysis: Optional[bool]) -> Dict[str, Any]:
    sh = qs.shared_paths(outdir)
    if run_processing is None:
        run_processing = bool(config.get("run_processing", True)) if config else None
    if run_analysis is None:
        run_analysis = bool(config.get("run_analysis", True)) if config else None
    cpath = counts_path if counts_path and os.path.exists(counts_path) else sh["deseq2"]
    counts, vst = None, None
    if os.path.exists(cpath):
        try:
            counts = qs.read_matrix(cpath)
            counts = counts[[c for c in counts.columns if c in set(sheet["Name"])]]
            if counts.shape[1] == 0:
                counts = None
        except Exception as exc:
            print(f"qc_report: could not read count matrix {cpath}: {exc}", file=sys.stderr)
    vpath = os.path.join(outdir, "analysis", "deseq2", "vst.tsv")
    if os.path.exists(vpath):
        try:
            vst = qs.read_matrix(vpath)
            vst = vst[[c for c in vst.columns if c in set(sheet["Name"])]]
        except Exception as exc:
            print(f"qc_report: could not read {vpath}: {exc}", file=sys.stderr)
    dpath = os.path.join(outdir, "analysis", "deseq2", "deg_summary.tsv")
    deg = _read_tsv(dpath)
    if deg is not None and not {"comparison", "n_up", "n_down"} <= set(deg.columns):
        print(f"qc_report: {dpath} lacks comparison/n_up/n_down columns; ignored", file=sys.stderr)
        deg = None
    if deg is not None and deg.empty:
        deg = None
    conds = list(dict.fromkeys(sheet["Condition"]))
    order = list(config.get("condition_order", []) or [])
    if order and set(order) == set(conds):
        conds = order
    bws = sorted(glob.glob(os.path.join(outdir, "tracks", "**", "*.bw"), recursive=True))
    return {"outdir": outdir, "sheet": sheet, "config": config or {}, "conditions": conds, "counts": counts,
            "counts_path": cpath, "vst": vst, "vst_path": vpath, "deg": deg, "deg_path": dpath,
            "scale": _read_tsv(os.path.join(outdir, "tracks", "scale_factors.tsv")), "bigwigs": bws,
            "run_processing": run_processing, "run_analysis": run_analysis}


def build_report(qc: pd.DataFrame, ctx: Dict[str, Any], thresholds, title: str = "TTpipe QC report",
                 banner: str = "") -> str:
    colours = cond_colours(ctx["conditions"])
    skips: Skips = []
    statuses = {r["sample"]: classify_sample({k: (None if not isnum(r.get(k)) else float(r[k])) for k in thresholds}, thresholds)
                for _, r in qc.iterrows()}
    overall = worst_status(st["overall"] for st in statuses.values())
    sections = [render_header(ctx, qc), render_summary(qc, thresholds, statuses, skips, ctx)]
    sections.append(render_alignment(qc, ctx, thresholds, skips))
    sections.append(render_single_metric("dup", "Duplication", "Picard PERCENT_DUPLICATION (dashed/dotted lines: warn/fail).",
                                         "duplicate_pct", "duplicates (%)", qc, ctx, thresholds, colours, skips))
    sections.append(render_strand(qc, ctx, thresholds, colours, skips))
    sections.append(render_single_metric(
        "rrna", "rRNA content", "% of counted reads in genes annotated rRNA / Mt_rRNA. mm10 references hold few rRNA loci, "
        "so this underestimates real rRNA contamination.", "rrna_pct", "rRNA reads (%)", qc, ctx, thresholds, colours, skips))
    sections.append(render_correlation(qc, ctx, thresholds, colours, skips))
    sections.append(render_pca(ctx, colours, skips))
    sections.append(render_deg(ctx, skips))
    sections.append(render_tracks(ctx, skips))
    sections.append(render_not_assessed())
    sections.append(render_methods(thresholds))
    # skipped-sections box goes right after the run header
    if skips:
        li = "".join(f"<li><b>{esc(a)}</b>: {esc(b)}</li>" for a, b in skips)
        sk = f'<section id="skipped"><h2>Skipped or partial sections</h2><div class="note"><ul>{li}</ul></div></section>'
    else:
        sk = '<section id="skipped"><h2>Skipped or partial sections</h2><p class="note">None: every section had data.</p></section>'
    sections.insert(1, sk)
    nav = [("run", "Run"), ("skipped", "Skipped"), ("summary", "Per-sample QC"), ("alignment", "Alignment"), ("dup", "Duplication"),
           ("strand", "Strandedness"), ("rrna", "rRNA"), ("correlation", "Correlation"), ("pca", "PCA"),
           ("deg", "DEGs"), ("tracks", "Tracks"), ("notassessed", "Not assessed"), ("methods", "Methods")]
    navhtml = "".join(f'<a href="#{i}">{esc(t)}</a>' for i, t in nav)
    head_badge = _badge(overall, {"pass": "All samples PASS", "warn": "Some samples WARN", "fail": "Some samples FAIL",
                                  "na": "No QC data"}[overall], "overall")
    now = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    ban = f'<div class="banner" role="note">{esc(banner)}</div>' if banner else ""
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light dark"><title>{esc(title)}</title>
<style>{CSS}</style>
<script>{get_plotlyjs()}</script></head>
<body><div class="wrap">
<header><div><h1>{esc(title)}</h1><div class="sub">{len(qc)} sample(s) &middot; generated {now}</div></div>
<div role="status" aria-label="Overall QC status">{head_badge}</div></header>
{ban}
<nav aria-label="Sections">{navhtml}</nav>
{''.join(sections)}
</div><script>{JS}</script></body></html>
"""


def _flag(v: Optional[str]) -> Optional[bool]:
    return None if v is None else str(v).lower() in ("1", "true", "yes")


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Generate the static TTpipe QC report.")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--config", nargs="*", default=[])
    ap.add_argument("--sample-qc", default=None, help="sample_qc.tsv (default OUTDIR/qc/sample_qc.tsv)")
    ap.add_argument("--counts", default=None, help="count matrix (analysis_counts); default OUTDIR/counts/counts_deseq2.txt")
    ap.add_argument("--run-processing", default=None, help="1/0; default from the config")
    ap.add_argument("--run-analysis", default=None, help="1/0; default from the config")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="TTpipe QC report")
    ap.add_argument("--banner", default="", help="text shown in a banner at the top (tests use it to label synthetic data)")
    a = ap.parse_args(argv)
    config = qs.load_config(a.config)
    sheet = qs.read_samplesheet(a.samplesheet)
    ctx = load_context(a.outdir, sheet, config, a.counts, _flag(a.run_processing), _flag(a.run_analysis))
    qpath = a.sample_qc or os.path.join(a.outdir, "qc", "sample_qc.tsv")
    qc = _read_tsv(qpath)
    if qc is None or "sample" not in qc.columns:
        print(f"qc_report: warning: {qpath} missing; computing the sample table directly", file=sys.stderr)
        qc = qs.build_sample_qc(a.outdir, sheet, config, a.counts)
    for c in qs.METRIC_COLUMNS:
        if c not in qc.columns:
            qc[c] = np.nan
    qc["sample"] = qc["sample"].astype(str)
    qc["condition"] = qc["condition"].astype(str)
    doc = build_report(qc, ctx, load_thresholds(config), a.title, a.banner)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(doc)
    print(f"qc_report: wrote {a.out} ({len(qc)} sample(s))", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
