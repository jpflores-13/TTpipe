#!/usr/bin/env python3
"""Static, self-contained TTpipe results explorer (one HTML file, plotly.js inlined; works from file://).

Usage:
    python scripts/results_explorer.py --degs OUT/analysis/deseq2/degs_*.csv \
        --gene-info OUT/counts/gene_info.tsv --norm OUT/analysis/deseq2/normalized_counts.tsv \
        --samples OUT/counts/samples_deseq2.txt --deg-summary OUT/analysis/deseq2/deg_summary.tsv \
        [--gsea Hallmarks=2H_vs_DMSO=OUT/analysis/gsea/Hallmarks_2H_vs_DMSO.tsv ...] \
        --alpha 0.01 --lfc 0 --out OUT/analysis/TTpipe_results_explorer.html

All data are embedded as compact JSON (arrays, not repeated keys); the controls are plain JavaScript, so the
file needs no network, no server and no framework. The module is split in two so the data preparation can be
tested without a browser:
    prepare_data(...)  tables -> a JSON-serialisable dict (no HTML)
    build_html(data)   dict -> the HTML document
Nothing is invented: a table that is missing (normalized counts, GSEA, gene annotation) is reported as
unavailable in the page, and the GSEA tab is absent when there are no GSEA results.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import os
import re
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

SIZE_WARN_MB = 40.0
PRIMARY_BIOTYPES = ["protein_coding", "lncRNA"]


# --------------------------------------------------------------------------
# Data preparation (no HTML here)
# --------------------------------------------------------------------------

def sig_round(v: Any, digits: int = 5) -> Optional[float]:
    """Round to `digits` significant digits; NaN / inf / None -> None (JSON null)."""
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    if math.isnan(f) or math.isinf(f):
        return None
    if f == 0.0:
        return 0.0
    return float(f"{f:.{digits - 1}e}")


def _col(values, digits: int) -> List[Optional[float]]:
    return [sig_round(v, digits) for v in values]


def contrast_name_from_path(path: str) -> str:
    """degs_<name>.csv -> <name> (the DEGS template in workflows/rules/common.smk)."""
    base = os.path.basename(path)
    m = re.fullmatch(r"degs_(.+)\.csv", base)
    return m.group(1) if m else os.path.splitext(base)[0]


def read_deg_table(path: str) -> pd.DataFrame:
    """One DESeq2 results table (gene_id first column; baseMean, log2FoldChange, padj used)."""
    df = pd.read_csv(path)
    if df.shape[1] == 0:
        raise ValueError(f"{path}: empty table")
    df = df.rename(columns={df.columns[0]: "gene_id"})
    df["gene_id"] = df["gene_id"].astype(str)
    for c in ("baseMean", "log2FoldChange", "padj"):
        if c not in df.columns:
            df[c] = np.nan
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df[["gene_id", "baseMean", "log2FoldChange", "padj"]]


def read_gene_info(path: Optional[str]) -> Optional[pd.DataFrame]:
    if not path or not os.path.exists(path):
        return None
    gi = pd.read_csv(path, sep="\t", dtype=str)
    if "gene_id" not in gi.columns:
        return None
    for c in ("gene_name", "gene_type"):
        if c not in gi.columns:
            gi[c] = ""
    gi = gi.drop_duplicates("gene_id").set_index("gene_id")
    return gi[["gene_name", "gene_type"]]


def read_samples(path: Optional[str]) -> Optional[pd.DataFrame]:
    if not path or not os.path.exists(path):
        return None
    s = pd.read_csv(path, sep="\t", dtype=str)
    if not {"sample", "condition"} <= set(s.columns):
        return None
    return s


def order_biotypes(types: List[str]) -> List[str]:
    """protein_coding, lncRNA first (when present), then the others by decreasing gene count, then name."""
    counts = pd.Series(types).value_counts()
    first = [t for t in PRIMARY_BIOTYPES if t in counts.index]
    rest = sorted((t for t in counts.index if t not in first), key=lambda t: (-int(counts[t]), t))
    return first + rest


def annotate_genes(gene_ids: List[str], gene_info: Optional[pd.DataFrame]) -> Tuple[List[str], List[str], int]:
    """(names, types, n_unannotated). gene_name falls back to gene_id, a missing gene_type to 'unknown'."""
    names, types, missing = [], [], 0
    for g in gene_ids:
        nm, ty = "", ""
        if gene_info is not None and g in gene_info.index:
            nm, ty = gene_info.at[g, "gene_name"], gene_info.at[g, "gene_type"]
            nm = "" if pd.isna(nm) else str(nm).strip()
            ty = "" if pd.isna(ty) else str(ty).strip()
        else:
            missing += 1
        names.append(nm if nm and nm != "NA" else g)
        types.append(ty if ty and ty != "NA" else "unknown")
    return names, types, missing


def read_gsea(path: str) -> List[list]:
    """[pathway, size, NES, pval, padj, leadingEdge] rows sorted by padj (NA last)."""
    df = pd.read_csv(path, sep="\t", dtype={"pathway": str, "leadingEdge": str})
    need = ["pathway", "size", "NES", "pval", "padj", "leadingEdge"]
    for c in need:
        if c not in df.columns:
            df[c] = np.nan
    for c in ("size", "NES", "pval", "padj"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.sort_values("padj", na_position="last", kind="stable")
    rows = []
    for p, sz, nes, pv, pa, le in zip(df["pathway"], df["size"], df["NES"], df["pval"], df["padj"], df["leadingEdge"]):
        rows.append([str(p), None if pd.isna(sz) else int(sz), sig_round(nes, 4), sig_round(pv, 4), sig_round(pa, 4),
                     "" if pd.isna(le) else str(le)])
    return rows


def prepare_data(degs: Dict[str, str], gene_info: Optional[str] = None, norm: Optional[str] = None,
                 samples: Optional[str] = None, gsea: Optional[Dict[str, Dict[str, str]]] = None,
                 deg_summary: Optional[str] = None, alpha: float = 0.01, lfc: float = 0.0,
                 condition_order: Optional[List[str]] = None, banner: str = "",
                 title: str = "TTpipe results explorer", generated: Optional[str] = None) -> Dict[str, Any]:
    """Read the pipeline tables and return the dict that is embedded in the page as JSON.

    degs: contrast name -> DESeq2 csv (insertion order = dropdown order). Genes whose baseMean is NA in a
    contrast are dropped from that contrast. gsea: collection -> contrast -> fgsea tsv (missing files skipped).
    """
    if not degs:
        raise ValueError("no DESeq2 result tables given")
    unavailable: List[str] = []
    tables = {name: read_deg_table(p) for name, p in degs.items()}
    for name, df in tables.items():
        tables[name] = df[df["baseMean"].notna()].drop_duplicates("gene_id")

    gene_ids: List[str] = list(dict.fromkeys(g for df in tables.values() for g in df["gene_id"]))
    index = {g: i for i, g in enumerate(gene_ids)}
    gi = read_gene_info(gene_info)
    if gi is None:
        unavailable.append("Gene annotation (gene_info.tsv) was not found: gene names equal gene IDs and the "
                           "biotype is 'unknown'.")
    names, types, n_unann = annotate_genes(gene_ids, gi)
    biotypes = order_biotypes(types) if types else []
    bt_index = {t: i for i, t in enumerate(biotypes)}

    contrasts = []
    for name, df in tables.items():
        contrasts.append({
            "name": name,
            "g": [index[g] for g in df["gene_id"]],
            "bm": _col(df["baseMean"], 5),
            "lfc": _col(df["log2FoldChange"], 5),
            "padj": _col(df["padj"], 4),
        })

    # sample table
    st = read_samples(samples)
    norm_block: Optional[Dict[str, Any]] = None
    sample_names: List[str] = []
    cond_of: Dict[str, str] = {}
    if st is not None:
        sample_names = list(st["sample"])
        cond_of = dict(zip(st["sample"], st["condition"]))
    else:
        unavailable.append("Sample table (samples_deseq2.txt) was not found: samples are not grouped by condition.")
    if norm and os.path.exists(norm):
        nm = pd.read_csv(norm, sep="\t", index_col=0)
        nm.index = nm.index.astype(str)
        nm = nm[~nm.index.duplicated()]
        nm.columns = [str(c) for c in nm.columns]
        cols = [s for s in sample_names if s in nm.columns] if sample_names else list(nm.columns)
        if not cols:
            cols = list(nm.columns)
        sub = nm.reindex(index=gene_ids, columns=cols)
        rows = []
        for vals in sub.to_numpy(dtype=float):
            rows.append(None if np.isnan(vals).all() else [sig_round(v, 4) for v in vals])
        norm_block = {"samples": cols, "cond": [cond_of.get(s, "n/a") for s in cols], "v": rows}
        if not sample_names:
            sample_names = cols
    else:
        unavailable.append("Normalized counts (normalized_counts.tsv) were not found: the per-gene count plot is "
                           "unavailable.")

    conditions = list(dict.fromkeys(cond_of.get(s, "n/a") for s in sample_names)) if sample_names else []
    if condition_order:
        known = [c for c in condition_order if c in conditions]
        conditions = known + [c for c in conditions if c not in known]

    gsea_block: Dict[str, Dict[str, list]] = {}
    for coll, per in (gsea or {}).items():
        for cname, path in per.items():
            if path and os.path.exists(path):
                rows = read_gsea(path)
                if rows:
                    gsea_block.setdefault(coll, {})[cname] = rows

    summary_rows = []
    if deg_summary and os.path.exists(deg_summary):
        sm = pd.read_csv(deg_summary, sep="\t")
        if {"comparison", "n_up", "n_down", "n_tested"} <= set(sm.columns):
            summary_rows = [[str(r.comparison), int(r.n_up), int(r.n_down), int(r.n_tested)] for r in sm.itertuples()]

    counts = pd.Series(types).value_counts() if types else pd.Series(dtype=int)
    return {
        "meta": {
            "title": title, "banner": banner or "",
            "generated": generated or _dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "alpha": float(alpha), "lfc": float(lfc),
            "samples": sample_names, "conditions": conditions,
            "n_genes": len(gene_ids), "n_unannotated": n_unann,
            "unavailable": unavailable, "deg_summary": summary_rows,
        },
        "biotypes": biotypes,
        "biotype_counts": [int(counts.get(t, 0)) for t in biotypes],
        "genes": {"id": gene_ids, "name": names, "type": [bt_index[t] for t in types]},
        "contrasts": contrasts,
        "norm": norm_block,
        "gsea": gsea_block,
    }


def to_json(data: Dict[str, Any]) -> str:
    """Compact JSON safe to place inside a <script type="application/json"> element."""
    s = json.dumps(data, separators=(",", ":"), allow_nan=False, ensure_ascii=False)
    return s.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")


# --------------------------------------------------------------------------
# HTML
# --------------------------------------------------------------------------

CSS = """
:root{--bg:#fff;--fg:#1b1f24;--muted:#5b6570;--card:#f6f8fa;--border:#d0d7de;--accent:#0072B2;--warn:#8a5a00;--warnbg:#fbeccb;--na:#5b6570}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#0f1318;--fg:#e6e9ed;--muted:#9aa5b1;--card:#171c23;
--border:#2f3741;--accent:#56B4E9;--warn:#f4c04d;--warnbg:#4a3a0a;--na:#9aa5b1}}
:root[data-theme="dark"]{--bg:#0f1318;--fg:#e6e9ed;--muted:#9aa5b1;--card:#171c23;--border:#2f3741;--accent:#56B4E9;--warn:#f4c04d;--warnbg:#4a3a0a;--na:#9aa5b1}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}
.wrap{max-width:1250px;margin:0 auto;padding:0 16px 64px}
header{padding:24px 0 8px}h1{font-size:1.6rem;margin:0}h2{font-size:1.2rem;margin:28px 0 8px;border-bottom:1px solid var(--border);padding-bottom:6px}
.sub{color:var(--muted);font-size:.9rem}
.note{background:var(--card);border-left:4px solid var(--na);padding:8px 12px;border-radius:4px;color:var(--muted);margin:8px 0}
.banner{background:var(--warnbg);color:var(--warn);border:1px solid var(--warn);border-radius:6px;padding:8px 12px;margin:12px 0;font-weight:600}
.err{background:var(--warnbg);color:var(--warn);border:1px solid var(--warn);border-radius:6px;padding:8px 12px;margin:12px 0}
details{margin:8px 0}summary{cursor:pointer;font-weight:600}
.tabs{display:flex;gap:4px;border-bottom:1px solid var(--border);margin-top:16px}
.tabs button{all:unset;cursor:pointer;padding:8px 16px;border:1px solid transparent;border-bottom:none;border-radius:6px 6px 0 0;font-weight:600;color:var(--muted)}
.tabs button[aria-selected="true"]{color:var(--fg);border-color:var(--border);background:var(--card)}
.tabs button:focus-visible{outline:2px solid var(--accent)}
.controls{display:flex;flex-wrap:wrap;gap:12px 24px;align-items:flex-end;background:var(--card);border:1px solid var(--border);border-radius:8px;padding:10px 14px;margin:12px 0}
.ctl{display:flex;flex-direction:column;gap:2px;font-size:.85rem;color:var(--muted)}
.ctl .row{display:flex;gap:8px;align-items:center;color:var(--fg)}
select,input[type=text],input[type=number]{background:var(--bg);color:var(--fg);border:1px solid var(--border);border-radius:4px;padding:4px 6px;font:inherit}
input[type=number]{width:90px}input[type=range]{width:140px}
button.btn{background:var(--accent);color:#fff;border:none;border-radius:4px;padding:5px 12px;font:inherit;cursor:pointer}
.bio{display:flex;flex-wrap:wrap;gap:4px 16px;align-items:center}.bio label{white-space:nowrap}
.bio-others{display:flex;flex-wrap:wrap;gap:2px 14px;margin-top:6px;max-height:160px;overflow:auto}
.grid{display:grid;grid-template-columns:3fr 2fr;gap:12px}
@media (max-width:900px){.grid{grid-template-columns:1fr}}
#plot,#countplot{min-height:480px}#countplot{min-height:360px}
.stat{font-weight:600;margin:6px 0}
.up{color:#E72A8A}.down{color:#66A61E}
.tablewrap{overflow:auto;max-height:520px;border:1px solid var(--border);border-radius:8px}
table{border-collapse:collapse;width:100%;font-size:.88rem}
th,td{padding:5px 10px;text-align:right;border-bottom:1px solid var(--border);white-space:nowrap}
td.l,th.l{text-align:left}
thead th{background:var(--card);position:sticky;top:0}
th button{all:unset;cursor:pointer;font-weight:600}th button:focus-visible{outline:2px solid var(--accent)}
th[aria-sort="ascending"] button::after{content:" \\25B2";font-size:.7em}
th[aria-sort="descending"] button::after{content:" \\25BC";font-size:.7em}
td.le{max-width:420px;overflow:hidden;text-overflow:ellipsis;cursor:pointer;text-align:left}
td.le.open{white-space:normal;word-break:break-all}
code{background:var(--card);padding:1px 5px;border-radius:4px}
dl{display:grid;grid-template-columns:max-content 1fr;gap:4px 16px}dt{font-weight:600}dd{margin:0}
.js-plotly-plot .plotly text{fill:var(--fg)}
[hidden]{display:none!important}
"""

JS = r"""
(function(){
'use strict';
var D = JSON.parse(document.getElementById('ttdata').textContent);
var $ = function(id){ return document.getElementById(id); };
var COL = {up:'#E72A8A', down:'#66A61E', ns:'#9aa0a6'};
var PAL = ['#0072B2','#E69F00','#009E73','#CC79A7','#56B4E9','#D55E00','#F0E442'];
var NG = D.genes.id.length;
var S = {ci:0, plot:'volcano', alpha:D.meta.alpha, lfc:D.meta.lfc, on:[], byBio:false, topOn:false, topN:10,
         sel:null, cat:'sig', text:'', sortKey:'padj', sortAsc:true, shown:200, view:null, clickBound:false,
         gColl:null, gCon:null, gText:'', gShown:300};
D.biotypes.forEach(function(){ S.on.push(true); });
var lowName = null, lowId = null;

function esc(s){ return String(s).replace(/[&<>"']/g, function(c){ return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]; }); }
function fmtN(v){ return v===null ? 'NA' : (Math.abs(v) >= 1000 || (v !== 0 && Math.abs(v) < 0.001) ? v.toExponential(3) : String(+v.toPrecision(4))); }
function fmtP(v){ return v===null ? 'NA' : v.toExponential(2); }
function fg(){ try { return getComputedStyle(document.body).color; } catch(e){ return '#222'; } }
function gridCol(){ return 'rgba(128,128,128,0.25)'; }
function catName(m){ return m===1 ? 'Up' : (m===2 ? 'Down' : 'n.s.'); }

function cls(c, k){
  var p = c.padj[k], l = c.lfc[k];
  if(p === null || l === null || p >= S.alpha || Math.abs(l) <= S.lfc) return 0;
  return l > 0 ? 1 : 2;
}
function compute(){
  var c = D.contrasts[S.ci], n = c.g.length, cat = new Int8Array(n), keep = new Uint8Array(n);
  var up = 0, down = 0, ns = 0, hidden = 0, naP = 0;
  for(var k = 0; k < n; k++){
    var gi = c.g[k];
    if(!S.on[D.genes.type[gi]]){ hidden++; continue; }
    keep[k] = 1;
    var m = cls(c, k); cat[k] = m;
    if(m === 1) up++; else if(m === 2) down++; else ns++;
    if(c.padj[k] === null) naP++;
  }
  S.view = {c:c, cat:cat, keep:keep, up:up, down:down, ns:ns, hidden:hidden, naP:naP, pos:null};
  return S.view;
}
function posOf(c){
  if(!c._pos){ c._pos = {}; for(var k = 0; k < c.g.length; k++) c._pos[c.g[k]] = k; }
  return c._pos;
}
function xy(c, k){
  var l = c.lfc[k];
  if(S.plot === 'volcano'){
    if(l === null) return null;
    var p = c.padj[k];
    return [l, p === null ? 0 : -Math.log10(Math.max(p, 1e-300))];
  }
  if(l === null || !(c.bm[k] > 0)) return null;
  return [Math.log10(c.bm[k]), l];
}
function hoverText(c, k){
  var gi = c.g[k];
  return '<b>' + esc(D.genes.name[gi]) + '</b><br>' + esc(D.genes.id[gi]) + '<br>' + esc(D.biotypes[D.genes.type[gi]]) +
    '<br>baseMean ' + fmtN(c.bm[k]) + '<br>log2FC ' + fmtN(c.lfc[k]) + '<br>padj ' + fmtP(c.padj[k]);
}
function bioGroup(t){ return t < PAL.length ? t : PAL.length; }   // first 7 biotypes get a colour, the rest are grouped

function drawPlot(){
  var v = compute(), c = v.c, groups = [], gmap = {};
  function grp(key, name, color){
    if(!gmap[key]){ gmap[key] = {name:name, color:color, x:[], y:[], text:[], cd:[]}; groups.push(gmap[key]); }
    return gmap[key];
  }
  var order = S.byBio ? null : [['ns','n.s.',COL.ns],['down','Down',COL.down],['up','Up',COL.up]];
  if(order){ order.forEach(function(o){ grp(o[0], o[1], o[2]); }); }
  else {
    D.biotypes.forEach(function(b, i){
      if(i < PAL.length) grp('b' + i, b, PAL[i]); else if(i === PAL.length) grp('b' + i, 'other biotypes', '#9aa0a6');
    });
  }
  var shown = 0;
  for(var k = 0; k < c.g.length; k++){
    if(!v.keep[k]) continue;
    var p = xy(c, k); if(p === null) continue;
    var g;
    if(S.byBio){ var bi = bioGroup(D.genes.type[c.g[k]]); g = gmap['b' + bi]; if(!g) continue; }
    else g = gmap[v.cat[k] === 1 ? 'up' : (v.cat[k] === 2 ? 'down' : 'ns')];
    g.x.push(p[0]); g.y.push(p[1]); g.text.push(hoverText(c, k)); g.cd.push(c.g[k]); shown++;
  }
  var traces = groups.filter(function(g){ return g.x.length > 0; }).map(function(g){
    return {type:'scatter', mode:'markers', name:g.name + ' (' + g.x.length + ')', x:g.x, y:g.y, text:g.text, customdata:g.cd,
            hoverinfo:'text', marker:{size:S.byBio ? 5 : 4, color:g.color, opacity:g.name === 'n.s.' ? 0.45 : 0.75}};
  });
  if(S.topOn && S.topN > 0){
    var cand = [];
    for(var j = 0; j < c.g.length; j++){ if(v.keep[j] && v.cat[j] !== 0 && xy(c, j) !== null) cand.push(j); }
    cand.sort(function(a, b){ return (c.padj[a] - c.padj[b]) || (Math.abs(c.lfc[b]) - Math.abs(c.lfc[a])); });
    cand = cand.slice(0, S.topN);
    traces.push({type:'scatter', mode:'text', name:'top ' + cand.length, showlegend:false, hoverinfo:'skip',
      x:cand.map(function(j){ return xy(c, j)[0]; }), y:cand.map(function(j){ return xy(c, j)[1]; }),
      text:cand.map(function(j){ return D.genes.name[c.g[j]]; }), textposition:'top center', textfont:{size:10, color:fg()}});
  }
  var selMsg = '';
  if(S.sel !== null){
    var kk = posOf(c)[S.sel], pp = (kk === undefined || !v.keep[kk]) ? null : xy(c, kk);
    if(pp){
      traces.push({type:'scatter', mode:'markers+text', name:D.genes.name[S.sel], showlegend:false, x:[pp[0]], y:[pp[1]],
        text:[D.genes.name[S.sel]], textposition:'top center', textfont:{size:12, color:fg()}, hoverinfo:'text',
        hovertext:[hoverText(c, kk)], customdata:[S.sel],
        marker:{size:12, color:'rgba(0,0,0,0)', line:{width:2.5, color:fg()}}});
    } else selMsg = ' Selected gene ' + D.genes.name[S.sel] + ' is not shown in this contrast (absent, all-NA, or hidden by the biotype filter).';
  }
  var shapes = [], vol = S.plot === 'volcano';
  if(vol){
    shapes.push({type:'line', xref:'paper', x0:0, x1:1, yref:'y', y0:-Math.log10(S.alpha), y1:-Math.log10(S.alpha), line:{dash:'dot', width:1, color:'rgba(128,128,128,.7)'}});
    if(S.lfc > 0) [-S.lfc, S.lfc].forEach(function(x){ shapes.push({type:'line', yref:'paper', y0:0, y1:1, xref:'x', x0:x, x1:x, line:{dash:'dot', width:1, color:'rgba(128,128,128,.7)'}}); });
  } else if(S.lfc > 0){
    [-S.lfc, S.lfc].forEach(function(y){ shapes.push({type:'line', xref:'paper', x0:0, x1:1, yref:'y', y0:y, y1:y, line:{dash:'dot', width:1, color:'rgba(128,128,128,.7)'}}); });
  }
  var layout = {
    title:{text:esc(c.name) + (vol ? ' (volcano)' : ' (MA)'), font:{size:14}},
    xaxis:{title:{text:vol ? 'log2 fold change' : 'log10(baseMean)'}, gridcolor:gridCol(), zeroline:false},
    yaxis:{title:{text:vol ? '-log10(padj)' : 'log2 fold change'}, gridcolor:gridCol(), zeroline:!vol},
    shapes:shapes, hovermode:'closest', paper_bgcolor:'rgba(0,0,0,0)', plot_bgcolor:'rgba(0,0,0,0)',
    font:{color:fg()}, margin:{l:60, r:10, t:40, b:50}, legend:{orientation:'h', y:-0.18}, height:520
  };
  $('counts').innerHTML = '<span class="up">Up: ' + v.up + '</span> &middot; <span class="down">Down: ' + v.down +
    '</span> &middot; n.s.: ' + v.ns + ' &nbsp;(' + (v.up + v.down + v.ns) + ' genes shown, ' + v.hidden + ' hidden by biotype filter; ' +
    v.naP + ' with NA padj counted as n.s.' + (vol ? ', drawn at y = 0' : '') + ')' + esc(selMsg);
  var pr = Plotly.react('plot', traces, layout, {responsive:true, displaylogo:false});
  return pr.then(function(){
    var gd = $('plot');
    if(!S.clickBound){
      S.clickBound = true;
      gd.on('plotly_click', function(ev){
        var pt = ev && ev.points && ev.points[0];
        if(pt && pt.customdata !== undefined && pt.customdata !== null) selectGene(pt.customdata);
      });
    }
    drawTable();
  });
}

function drawCounts(){
  var note = $('countnote'), box = $('countplot');
  function say(t){ note.textContent = t; note.hidden = false; box.hidden = true; return Promise.resolve(); }
  if(!D.norm) return say('Normalized counts are unavailable in this run (normalized_counts.tsv was not found).');
  if(S.sel === null) return say('Click a point in the plot, or search for a gene, to see its normalized counts per sample.');
  var row = D.norm.v[S.sel];
  if(!row) return say('No normalized counts for ' + D.genes.name[S.sel] + ' (it is not in normalized_counts.tsv).');
  note.hidden = true; box.hidden = false;
  var conds = D.meta.conditions.length ? D.meta.conditions : Array.from(new Set(D.norm.cond));
  var traces = conds.map(function(cn, i){
    var xs = [], ys = [];
    D.norm.samples.forEach(function(s, j){ if(D.norm.cond[j] === cn){ xs.push(s); ys.push(row[j]); } });
    return {type:'bar', name:cn, x:xs, y:ys, marker:{color:PAL[i % PAL.length]}, hovertemplate:'%{x}<br>%{y}<extra>' + esc(cn) + '</extra>'};
  }).filter(function(t){ return t.x.length > 0; });
  var layout = {title:{text:esc(D.genes.name[S.sel]) + ' (' + esc(D.genes.id[S.sel]) + '): normalized counts', font:{size:13}},
    yaxis:{title:{text:'DESeq2 normalized count'}, gridcolor:gridCol()}, xaxis:{type:'category', gridcolor:gridCol()},
    paper_bgcolor:'rgba(0,0,0,0)', plot_bgcolor:'rgba(0,0,0,0)', font:{color:fg()}, margin:{l:60, r:10, t:40, b:90},
    legend:{orientation:'h', y:-0.35}, barmode:'group', height:360};
  return Plotly.react('countplot', traces, layout, {responsive:true, displaylogo:false});
}

function selectGene(gi){
  S.sel = gi;
  $('search').value = D.genes.name[gi];
  var a = drawPlot(), b = drawCounts();
  return Promise.all([a, b]);
}

/* ---- search ---- */
function buildLower(){
  if(lowName) return;
  lowName = D.genes.name.map(function(s){ return s.toLowerCase(); });
  lowId = D.genes.id.map(function(s){ return s.toLowerCase(); });
}
function matches(q, limit){
  buildLower();
  q = q.trim().toLowerCase();
  var out = [];
  if(!q) return out;
  for(var i = 0; i < NG && out.length < limit; i++){
    if(lowName[i].indexOf(q) === 0 || lowId[i].indexOf(q) === 0) out.push(i);
  }
  return out;
}
function onSearchInput(){
  var q = $('search').value, dl = $('hits'), m = matches(q, 30);
  dl.innerHTML = m.map(function(i){ return '<option value="' + esc(D.genes.name[i]) + '">' + esc(D.genes.id[i]) + '</option>'; }).join('');
}
function onSearchCommit(){
  var q = $('search').value.trim().toLowerCase(), st = $('searchmsg');
  if(!q){ S.sel = null; st.textContent = ''; drawPlot(); drawCounts(); return; }
  buildLower();
  var hit = -1;
  for(var i = 0; i < NG; i++){ if(lowName[i] === q || lowId[i] === q){ hit = i; break; } }
  if(hit < 0){ var m = matches(q, 1); if(m.length) hit = m[0]; }
  if(hit < 0){ st.textContent = 'No gene starts with "' + $('search').value + '".'; return; }
  st.textContent = '';
  selectGene(hit);
}

/* ---- results table ---- */
var COLS = [['name','Gene','l'],['id','gene_id','l'],['type','Biotype','l'],['bm','baseMean',''],['lfc','log2FC',''],['padj','padj',''],['cat','Call','l']];
function rowsFiltered(){
  var v = S.view, c = v.c, out = [], t = S.text.trim().toLowerCase();
  for(var k = 0; k < c.g.length; k++){
    if(!v.keep[k]) continue;
    var m = v.cat[k];
    if(S.cat === 'sig' && m === 0) continue;
    if(S.cat === 'up' && m !== 1) continue;
    if(S.cat === 'down' && m !== 2) continue;
    if(S.cat === 'ns' && m !== 0) continue;
    if(t){
      var gi = c.g[k];
      if((D.genes.name[gi] + ' ' + D.genes.id[gi] + ' ' + D.biotypes[D.genes.type[gi]]).toLowerCase().indexOf(t) < 0) continue;
    }
    out.push(k);
  }
  var key = S.sortKey, dir = S.sortAsc ? 1 : -1;
  function val(k){
    var gi = c.g[k];
    switch(key){
      case 'name': return D.genes.name[gi].toLowerCase();
      case 'id': return D.genes.id[gi].toLowerCase();
      case 'type': return D.biotypes[D.genes.type[gi]].toLowerCase();
      case 'bm': return c.bm[k];
      case 'lfc': return c.lfc[k];
      case 'padj': return c.padj[k];
      default: return v.cat[k];
    }
  }
  out.sort(function(a, b){
    var x = val(a), y = val(b);
    if(x === null && y === null) return 0;
    if(x === null) return 1;
    if(y === null) return -1;
    return (x < y ? -1 : (x > y ? 1 : 0)) * dir;
  });
  return out;
}
function drawTable(){
  var v = S.view, c = v.c, rows = rowsFiltered();
  S.rows = rows;
  var head = '<tr>' + COLS.map(function(col){
    var so = S.sortKey === col[0] ? (S.sortAsc ? 'ascending' : 'descending') : 'none';
    return '<th class="' + col[2] + '" aria-sort="' + so + '"><button type="button" data-sort="' + col[0] + '">' + col[1] + '</button></th>';
  }).join('') + '</tr>';
  var n = Math.min(rows.length, S.shown), body = [];
  for(var i = 0; i < n; i++){
    var k = rows[i], gi = c.g[k], m = v.cat[k];
    body.push('<tr data-gi="' + gi + '"><td class="l">' + esc(D.genes.name[gi]) + '</td><td class="l">' + esc(D.genes.id[gi]) +
      '</td><td class="l">' + esc(D.biotypes[D.genes.type[gi]]) + '</td><td>' + fmtN(c.bm[k]) + '</td><td>' + fmtN(c.lfc[k]) +
      '</td><td>' + fmtP(c.padj[k]) + '</td><td class="l ' + (m === 1 ? 'up' : (m === 2 ? 'down' : '')) + '">' + catName(m) + '</td></tr>');
  }
  $('tbl').innerHTML = '<thead>' + head + '</thead><tbody>' + body.join('') + '</tbody>';
  $('tblinfo').textContent = 'Showing ' + n + ' of ' + rows.length + ' genes (' + esc(c.name) + '; padj < ' + S.alpha + ', |log2FC| > ' + S.lfc + ').';
  $('more').hidden = rows.length <= n;
}
function currentCsv(){
  var c = S.view.c, v = S.view, rows = S.rows || rowsFiltered();
  function q(s){ s = String(s); return /[",\n\r]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; }
  var lines = ['gene_id,gene_name,gene_type,baseMean,log2FoldChange,padj,call,contrast'];
  rows.forEach(function(k){
    var gi = c.g[k];
    lines.push([q(D.genes.id[gi]), q(D.genes.name[gi]), q(D.biotypes[D.genes.type[gi]]), c.bm[k] === null ? 'NA' : c.bm[k],
      c.lfc[k] === null ? 'NA' : c.lfc[k], c.padj[k] === null ? 'NA' : c.padj[k], catName(v.cat[k]), q(c.name)].join(','));
  });
  return lines.join('\n') + '\n';
}
function downloadCsv(){
  var blob = new Blob([currentCsv()], {type:'text/csv;charset=utf-8'});
  var a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'TTpipe_' + D.contrasts[S.ci].name + '_shown.csv';
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  setTimeout(function(){ URL.revokeObjectURL(a.href); }, 1000);
}

/* ---- GSEA ---- */
function gseaRows(){
  var rows = (D.gsea[S.gColl] || {})[S.gCon] || [], t = S.gText.trim().toLowerCase();
  return t ? rows.filter(function(r){ return r[0].toLowerCase().indexOf(t) >= 0 || r[5].toLowerCase().indexOf(t) >= 0; }) : rows;
}
function drawGsea(){
  var rows = gseaRows(), n = Math.min(rows.length, S.gShown), body = [];
  for(var i = 0; i < n; i++){
    var r = rows[i], le = r[5], short = le.length > 60 ? le.slice(0, 60) + '\u2026' : le;
    body.push('<tr><td class="l">' + esc(r[0]) + '</td><td>' + (r[1] === null ? 'NA' : r[1]) + '</td><td>' + fmtN(r[2]) + '</td><td>' + fmtP(r[3]) +
      '</td><td>' + fmtP(r[4]) + '</td><td class="le" title="' + esc(le) + '" data-full="' + esc(le) + '" data-short="' + esc(short) + '">' + esc(short) + '</td></tr>');
  }
  $('gtbl').innerHTML = '<thead><tr><th class="l">pathway</th><th>size</th><th>NES</th><th>pval</th><th>padj</th><th class="l">leadingEdge (click to expand)</th></tr></thead><tbody>' +
    body.join('') + '</tbody>';
  $('ginfo').textContent = 'Showing ' + n + ' of ' + rows.length + ' pathways, sorted by padj (NA last).';
  $('gmore').hidden = rows.length <= n;
}

/* ---- biotype controls ---- */
function buildBiotypes(){
  var main = $('bio-main'), others = $('bio-others'), hm = '', ho = '';
  D.biotypes.forEach(function(b, i){
    var lab = '<label><input type="checkbox" data-bt="' + i + '" checked> ' + esc(b) + ' <span class="sub">(' + D.biotype_counts[i] + ')</span></label>';
    if(PRIMARY.indexOf(b) >= 0) hm += lab; else ho += lab;
  });
  main.innerHTML = hm;
  others.innerHTML = ho;
  $('bio-other-wrap').hidden = ho === '';
}
var PRIMARY = ['protein_coding', 'lncRNA'];
function syncGroupBox(){
  var tot = 0, on = 0;
  D.biotypes.forEach(function(b, i){ if(PRIMARY.indexOf(b) < 0){ tot++; if(S.on[i]) on++; } });
  var g = $('bio-group'); g.checked = tot > 0 && on === tot; g.indeterminate = on > 0 && on < tot;
}
function syncBoxes(){
  document.querySelectorAll('[data-bt]').forEach(function(cb){ cb.checked = !!S.on[+cb.getAttribute('data-bt')]; });
  syncGroupBox();
}

/* ---- slider / number linking ---- */
function bindCutoffs(){
  var an = $('alpha-n'), ar = $('alpha-r'), ln = $('lfc-n'), lr = $('lfc-r');
  function setAlpha(a){ if(!(a > 0) || a > 1) return; S.alpha = a; an.value = String(+a.toPrecision(3)); ar.value = Math.max(-10, Math.min(0, Math.log10(a))); S.shown = 200; drawPlot(); }
  function setLfc(l){ if(!(l >= 0)) return; S.lfc = l; ln.value = String(+l.toPrecision(3)); lr.value = Math.min(+lr.max, l); S.shown = 200; drawPlot(); }
  an.value = String(S.alpha); ar.value = Math.max(-10, Math.min(0, Math.log10(S.alpha))); ln.value = String(S.lfc); lr.value = Math.min(+lr.max, S.lfc);
  an.addEventListener('change', function(){ setAlpha(parseFloat(an.value)); });
  an.addEventListener('input', function(){ var a = parseFloat(an.value); if(a > 0 && a <= 1) setAlpha(a); });
  ar.addEventListener('input', function(){ setAlpha(Math.pow(10, parseFloat(ar.value))); });
  ln.addEventListener('change', function(){ setLfc(parseFloat(ln.value)); });
  ln.addEventListener('input', function(){ var l = parseFloat(ln.value); if(l >= 0) setLfc(l); });
  lr.addEventListener('input', function(){ setLfc(parseFloat(lr.value)); });
}

function init(){
  var sel = $('contrast');
  sel.innerHTML = D.contrasts.map(function(c, i){ return '<option value="' + i + '">' + esc(c.name) + '</option>'; }).join('');
  sel.addEventListener('change', function(){ S.ci = +sel.value; S.shown = 200; drawPlot(); drawCounts(); });
  document.querySelectorAll('input[name=ptype]').forEach(function(r){
    r.addEventListener('change', function(){ if(r.checked){ S.plot = r.value; drawPlot(); } });
  });
  bindCutoffs();
  buildBiotypes();
  $('bio-wrap').addEventListener('change', function(ev){
    var t = ev.target;
    if(t.id === 'bio-group'){
      D.biotypes.forEach(function(b, i){ if(PRIMARY.indexOf(b) < 0) S.on[i] = t.checked; });
      syncBoxes();
    } else if(t.hasAttribute('data-bt')){
      S.on[+t.getAttribute('data-bt')] = t.checked; syncGroupBox();
    } else return;
    S.shown = 200; drawPlot(); drawCounts();
  });
  $('bio-color').addEventListener('change', function(ev){ S.byBio = ev.target.checked; drawPlot(); });
  $('top-on').addEventListener('change', function(ev){ S.topOn = ev.target.checked; drawPlot(); });
  $('top-n').addEventListener('change', function(ev){ var n = parseInt(ev.target.value, 10); S.topN = n > 0 ? n : 0; drawPlot(); });
  $('search').addEventListener('input', onSearchInput);
  $('search').addEventListener('change', onSearchCommit);
  $('search').addEventListener('keydown', function(ev){ if(ev.key === 'Enter') onSearchCommit(); });
  $('tcat').addEventListener('change', function(ev){ S.cat = ev.target.value; S.shown = 200; drawTable(); });
  $('ttext').addEventListener('input', function(ev){ S.text = ev.target.value; S.shown = 200; drawTable(); });
  $('more').addEventListener('click', function(){ S.shown += 200; drawTable(); });
  $('csv').addEventListener('click', downloadCsv);
  $('tbl').addEventListener('click', function(ev){
    var b = ev.target.closest('button[data-sort]');
    if(b){
      var k = b.getAttribute('data-sort');
      if(S.sortKey === k) S.sortAsc = !S.sortAsc; else { S.sortKey = k; S.sortAsc = (k === 'name' || k === 'id' || k === 'type' || k === 'padj'); }
      drawTable(); return;
    }
    var tr = ev.target.closest('tr[data-gi]');
    if(tr) selectGene(+tr.getAttribute('data-gi'));
  });
  document.querySelectorAll('.tabs button').forEach(function(b){
    b.addEventListener('click', function(){
      document.querySelectorAll('.tabs button').forEach(function(o){ o.setAttribute('aria-selected', o === b ? 'true' : 'false'); });
      document.querySelectorAll('.panel').forEach(function(p){ p.hidden = p.id !== b.getAttribute('data-panel'); });
      if(b.getAttribute('data-panel') === 'p-results') window.dispatchEvent(new Event('resize'));
      if(b.getAttribute('data-panel') === 'p-gsea') drawGsea();
    });
  });
  var colls = Object.keys(D.gsea);
  if(colls.length){
    S.gColl = colls[0];
    var gc = $('gcoll'), gk = $('gcon');
    gc.innerHTML = colls.map(function(c){ return '<option>' + esc(c) + '</option>'; }).join('');
    function fillCon(){
      var cons = Object.keys(D.gsea[S.gColl]);
      gk.innerHTML = cons.map(function(c){ return '<option>' + esc(c) + '</option>'; }).join('');
      S.gCon = cons[0]; S.gShown = 300;
    }
    fillCon();
    gc.addEventListener('change', function(){ S.gColl = gc.value; fillCon(); drawGsea(); });
    gk.addEventListener('change', function(){ S.gCon = gk.value; S.gShown = 300; drawGsea(); });
    $('gtext').addEventListener('input', function(ev){ S.gText = ev.target.value; S.gShown = 300; drawGsea(); });
    $('gmore').addEventListener('click', function(){ S.gShown += 300; drawGsea(); });
    $('gtbl').addEventListener('click', function(ev){
      var td = ev.target.closest('td.le');
      if(td){ var open = td.classList.toggle('open'); td.textContent = open ? td.getAttribute('data-full') : td.getAttribute('data-short'); }
    });
  }
  return drawPlot().then(drawCounts);
}

/* ---- error reporting (visible, and in the DOM for automated checks) ---- */
function showErrors(){
  var errs = window.__ttErrors || [], box = $('errbox');
  if(errs.length){ box.hidden = false; box.textContent = 'JavaScript error(s): ' + errs.join(' | '); }
}

/* ---- optional self-test: open the file with #selftest to exercise the controls ---- */
function selfTest(){
  var out = [], fails = 0;
  function ok(name, cond, extra){ if(!cond) fails++; out.push((cond ? 'PASS ' : 'FAIL ') + name + (extra ? ' ' + extra : '')); }
  function fire(el, type){ el.dispatchEvent(new Event(type, {bubbles:true})); }
  function nPoints(){ return document.querySelectorAll('#plot .scatterlayer path.point').length; }
  var chain = Promise.resolve();
  function step(f){ chain = chain.then(f); }
  step(function(){
    ok('contrast options', $('contrast').options.length === D.contrasts.length, String($('contrast').options.length));
    ok('volcano points drawn', nPoints() > 0, String(nPoints()));
    ok('volcano axis', $('plot').layout.xaxis.title.text === 'log2 fold change');
    ok('counts shown', /Up: \d+/.test($('counts').textContent), $('counts').textContent.slice(0, 40));
    ok('table rows', document.querySelectorAll('#tbl tbody tr').length > 0, String(document.querySelectorAll('#tbl tbody tr').length));
  });
  step(function(){
    var r = document.querySelector('input[name=ptype][value=ma]'); r.checked = true; fire(r, 'change');
    return new Promise(function(res){ setTimeout(res, 50); });
  });
  step(function(){
    ok('MA axis', $('plot').layout.xaxis.title.text === 'log10(baseMean)');
    var before = S.view.up + S.view.down;
    $('alpha-n').value = '0.5'; fire($('alpha-n'), 'change');
    return new Promise(function(res){ setTimeout(function(){ ok('alpha 0.5 changes live counts', S.view.up + S.view.down >= before, before + ' -> ' + (S.view.up + S.view.down)); res(); }, 50); });
  });
  step(function(){
    $('alpha-n').value = '1e-12'; fire($('alpha-n'), 'change');
    return new Promise(function(res){ setTimeout(function(){ ok('alpha 1e-12 not more than 0.5 setting', S.view.up + S.view.down < 1e9); res(); }, 50); });
  });
  step(function(){
    var tot = S.view.up + S.view.down + S.view.ns;
    var cb = document.querySelector('[data-bt="0"]'); cb.checked = false; fire(cb, 'change');
    return new Promise(function(res){ setTimeout(function(){
      var tot2 = S.view.up + S.view.down + S.view.ns;
      ok('biotype filter hides genes', D.biotypes.length < 2 || tot2 < tot, tot + ' -> ' + tot2);
      cb.checked = true; fire(cb, 'change'); res(); }, 50); });
  });
  step(function(){
    var c = $('bio-color'); c.checked = true; fire(c, 'change');
    return new Promise(function(res){ setTimeout(function(){ ok('colour by biotype', nPoints() > 0 && $('plot').data.length >= 1, String($('plot').data.length) + ' traces'); c.checked = false; fire(c, 'change'); res(); }, 50); });
  });
  step(function(){
    var t = $('top-on'); t.checked = true; fire(t, 'change');
    return new Promise(function(res){ setTimeout(function(){
      var has = $('plot').data.some(function(d){ return d.mode === 'text'; });
      ok('top-N labels', has || S.view.up + S.view.down === 0); res(); }, 50); });
  });
  step(function(){
    var nm = D.genes.name[D.contrasts[S.ci].g[0]];
    $('search').value = nm.slice(0, Math.max(1, nm.length - 1)); fire($('search'), 'input');
    ok('search datalist', $('hits').options.length >= 1, String($('hits').options.length));
    $('search').value = nm; fire($('search'), 'change');
    return new Promise(function(res){ setTimeout(function(){
      ok('search selects gene', S.sel !== null && D.genes.name[S.sel] === nm, nm);
      ok('highlight trace', $('plot').data.some(function(d){ return d.mode === 'markers+text' && d.name === nm; }) || $('counts').textContent.indexOf('not shown') >= 0);
      if(D.norm) ok('count plot drawn', !$('countplot').hidden && document.querySelectorAll('#countplot .barlayer .point').length > 0, String(document.querySelectorAll('#countplot .barlayer .point').length) + ' bars');
      else ok('count plot unavailable message', !$('countnote').hidden && /unavailable/.test($('countnote').textContent));
      res(); }, 80); });
  });
  step(function(){
    var gi = D.contrasts[S.ci].g[Math.min(3, D.contrasts[S.ci].g.length - 1)];
    $('plot').emit('plotly_click', {points:[{customdata:gi}]});
    return new Promise(function(res){ setTimeout(function(){ ok('plotly_click selects gene', S.sel === gi); res(); }, 80); });
  });
  step(function(){
    var b = document.querySelector('#tbl th button[data-sort=lfc]'); b.click();
    var rows = S.rows, c = S.view.c, okSort = true, d = S.sortAsc ? 1 : -1;
    for(var i = 1; i < rows.length; i++){ if(c.lfc[rows[i]] !== null && c.lfc[rows[i - 1]] !== null && (c.lfc[rows[i]] - c.lfc[rows[i - 1]]) * d < 0) okSort = false; }
    ok('table sorts by log2FC (' + (S.sortAsc ? 'ascending' : 'descending') + ')', S.sortKey === 'lfc' && okSort, rows.length + ' rows');
    b.click(); d = S.sortAsc ? 1 : -1; okSort = true; rows = S.rows;
    for(var j = 1; j < rows.length; j++){ if((c.lfc[rows[j]] - c.lfc[rows[j - 1]]) * d < 0) okSort = false; }
    ok('second click reverses the sort', okSort);
    $('tcat').value = 'all'; fire($('tcat'), 'change');
    ok('table filter all', S.rows.length >= rows.length);
    $('ttext').value = 'zzzz_no_such_gene'; fire($('ttext'), 'input');
    ok('table text filter', S.rows.length === 0);
    $('ttext').value = ''; fire($('ttext'), 'input');
    var csv = currentCsv().trim().split('\n');
    ok('csv rows match table', csv.length - 1 === S.rows.length && csv[0].indexOf('gene_id,gene_name') === 0, (csv.length - 1) + ' rows');
  });
  if(D.contrasts.length > 1){
    step(function(){
      $('contrast').value = '1'; fire($('contrast'), 'change');
      return new Promise(function(res){ setTimeout(function(){ ok('second contrast drawn', S.ci === 1 && nPoints() > 0 && $('plot').layout.title.text.indexOf(D.contrasts[1].name) === 0); res(); }, 80); });
    });
  }
  if(Object.keys(D.gsea).length){
    step(function(){
      var tab = document.querySelector('.tabs button[data-panel=p-gsea]'); tab.click();
      ok('GSEA tab visible', !$('p-gsea').hidden && document.querySelectorAll('#gtbl tbody tr').length > 0, String(document.querySelectorAll('#gtbl tbody tr').length) + ' rows');
      var n0 = document.querySelectorAll('#gtbl tbody tr').length;
      $('gtext').value = 'zzzz_no_such_pathway'; fire($('gtext'), 'input');
      ok('GSEA text filter', document.querySelectorAll('#gtbl tbody tr').length === 0 && n0 > 0);
    });
  } else {
    step(function(){ ok('no GSEA tab when no GSEA data', !document.querySelector('.tabs button[data-panel=p-gsea]')); });
  }
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
  }).catch(function(e){ (window.__ttErrors = window.__ttErrors || []).push(String(e)); showErrors(); });
} catch(e){ (window.__ttErrors = window.__ttErrors || []).push(String(e)); showErrors(); }
window.addEventListener('error', showErrors);
})();
"""

BODY = """
<div class="wrap">
<header><h1>__TITLE__</h1><div class="sub">__SUB__</div></header>
__BANNER__
<div id="errbox" class="err" role="alert" hidden></div>
__UNAVAILABLE__
<details open><summary>How to read this page</summary>
<div class="note">
<p><b>Volcano</b>: x is log2 fold change (treatment vs reference of the chosen contrast), y is -log10 of the DESeq2 adjusted p-value (padj).
<b>MA</b>: x is log10(baseMean), y is log2 fold change. A gene is <span class="up">Up</span> when padj &lt; the padj cut-off and
log2FC &gt; the |log2FC| cut-off, <span class="down">Down</span> when padj is below the cut-off and log2FC is below minus the |log2FC| cut-off, and
n.s. otherwise. The two cut-offs start at the values of <code>alpha</code> and <code>lfc_threshold</code> from the config and can be changed here;
the pipeline files on disk (deg_summary.tsv, Sankey) always use the config values. If <code>lfc_threshold</code> &gt; 0 DESeq2 already tested
against that threshold, so padj and the |log2FC| filter here are not independent.</p>
<p>Genes with NA padj (DESeq2 independent filtering or outliers) are shown as n.s.; in the volcano they are drawn at y = 0. Genes with NA baseMean
are not in the file. Biotypes come from the annotation (<code>gene_type</code> in the GTF, GENCODE style): nothing is detected here, and a
lncRNA is a lncRNA only because the annotation says so. Click a point, a table row or use the search box to see the normalized counts of one gene.
There is no genome browser or track panel in this version.</p>
</div></details>
<div class="tabs" role="tablist">
<button type="button" role="tab" aria-selected="true" data-panel="p-results">Results</button>__GSEATAB__
</div>
<section id="p-results" class="panel">
<div class="controls">
<div class="ctl"><label for="contrast">Contrast</label><select id="contrast"></select></div>
<div class="ctl"><span>Plot</span><div class="row"><label><input type="radio" name="ptype" value="volcano" checked> Volcano</label>
<label><input type="radio" name="ptype" value="ma"> MA</label></div></div>
<div class="ctl"><label for="alpha-n">padj cut-off</label><div class="row"><input id="alpha-n" type="number" min="0" max="1" step="any">
<input id="alpha-r" type="range" min="-10" max="0" step="0.1" aria-label="padj cut-off (log10 scale)"></div></div>
<div class="ctl"><label for="lfc-n">|log2FC| cut-off</label><div class="row"><input id="lfc-n" type="number" min="0" step="any">
<input id="lfc-r" type="range" min="0" max="5" step="0.05" aria-label="log2FC cut-off"></div></div>
<div class="ctl"><label for="search">Find gene (symbol or gene_id, prefix)</label><div class="row"><input id="search" type="text" list="hits" autocomplete="off" placeholder="e.g. Gapdh">
<datalist id="hits"></datalist></div><span id="searchmsg" class="sub"></span></div>
<div class="ctl"><span>Labels and colours</span><div class="row"><label><input id="top-on" type="checkbox"> label top</label>
<input id="top-n" type="number" min="1" max="200" value="10" aria-label="number of genes to label"><span>genes</span>
<label><input id="bio-color" type="checkbox"> colour by biotype</label></div></div>
</div>
<div id="bio-wrap" class="controls"><div class="ctl"><span>Gene biotype (GENCODE gene_type)</span>
<div class="bio" id="bio-main"></div>
<div id="bio-other-wrap"><label><input type="checkbox" id="bio-group" checked> all other biotypes</label>
<details><summary>choose individual other biotypes</summary><div class="bio-others" id="bio-others"></div></details></div></div></div>
<div class="stat" id="counts"></div>
<div class="grid"><div id="plot"></div><div><div id="countnote" class="note"></div><div id="countplot" hidden></div></div></div>
<h2>Results table</h2>
<div class="controls">
<div class="ctl"><label for="tcat">Show</label><select id="tcat"><option value="sig" selected>Up and Down</option><option value="up">Up</option>
<option value="down">Down</option><option value="ns">n.s.</option><option value="all">all genes</option></select></div>
<div class="ctl"><label for="ttext">Filter (name, gene_id, biotype)</label><input id="ttext" type="text"></div>
<div class="ctl"><span>&nbsp;</span><button type="button" class="btn" id="csv">Download CSV of the rows matching the filters</button></div>
</div>
<div class="sub" id="tblinfo"></div>
<div class="tablewrap"><table id="tbl"></table></div>
<p><button type="button" class="btn" id="more" hidden>Show 200 more rows</button> <span class="sub">The CSV always contains every matching row, not just the rows displayed.</span></p>
</section>
__GSEAPANEL__
<pre id="selftest" hidden></pre>
</div>
"""

GSEA_PANEL = """<section id="p-gsea" class="panel" hidden>
<p class="note">fgsea preranked results written by the pipeline (genes ranked by the DESeq2 Wald statistic; <code>gsea_collections</code> in the config).
Sorted by padj. Rows with NA padj come last. The leadingEdge column is cut short; click a cell to expand it.</p>
<div class="controls">
<div class="ctl"><label for="gcoll">Collection</label><select id="gcoll"></select></div>
<div class="ctl"><label for="gcon">Contrast</label><select id="gcon"></select></div>
<div class="ctl"><label for="gtext">Filter (pathway or gene)</label><input id="gtext" type="text"></div>
</div>
<div class="sub" id="ginfo"></div>
<div class="tablewrap"><table id="gtbl"></table></div>
<p><button type="button" class="btn" id="gmore" hidden>Show 300 more rows</button></p>
</section>"""


def _esc(x: Any) -> str:
    import html
    return html.escape(str(x), quote=True)


def build_html(data: Dict[str, Any], plotly_js: Optional[str] = None) -> str:
    """The complete page. `plotly_js` defaults to the plotly.js bundled with the installed plotly package."""
    if plotly_js is None:
        from plotly.offline import get_plotlyjs
        plotly_js = get_plotlyjs()
    m = data["meta"]
    sub = (f"{len(m['samples'])} sample(s) &middot; {len(m['conditions'])} condition(s) ({_esc(', '.join(m['conditions'])) or 'n/a'}) &middot; "
           f"{len(data['contrasts'])} contrast(s) ({_esc(', '.join(c['name'] for c in data['contrasts']))}) &middot; "
           f"initial padj &lt; {_esc(m['alpha'])}, |log2FC| &gt; {_esc(m['lfc'])} &middot; {m['n_genes']} genes &middot; generated {_esc(m['generated'])}")
    banner = f'<div class="banner" role="note">{_esc(m["banner"])}</div>' if m.get("banner") else ""
    items = list(m.get("unavailable", []))
    if not data.get("gsea"):
        items.append("No GSEA results were found (gsea_collections is empty or the files are missing): there is no GSEA tab.")
    if m.get("n_unannotated"):
        items.append(f"{m['n_unannotated']} gene(s) are not in the gene table: their name is the gene_id and their biotype is 'unknown'.")
    unav = ""
    if items:
        unav = ('<details><summary>Not available in this file</summary><div class="note"><ul>' +
                "".join(f"<li>{_esc(i)}</li>" for i in items) + "</ul></div></details>")
    if m.get("deg_summary"):
        rows = "".join(f"<tr><td class='l'>{_esc(r[0])}</td><td>{r[1]}</td><td>{r[2]}</td><td>{r[3]}</td></tr>" for r in m["deg_summary"])
        unav += ("<details><summary>Pipeline DEG summary (config alpha and lfc_threshold, from deg_summary.tsv)</summary><div class='tablewrap'>"
                 "<table><thead><tr><th class='l'>comparison</th><th>n_up</th><th>n_down</th><th>n_tested</th></tr></thead>"
                 f"<tbody>{rows}</tbody></table></div></details>")
    gsea_tab = '<button type="button" role="tab" aria-selected="false" data-panel="p-gsea">GSEA</button>' if data.get("gsea") else ""
    gsea_panel = GSEA_PANEL if data.get("gsea") else ""
    body = (BODY.replace("__TITLE__", _esc(m["title"])).replace("__SUB__", sub).replace("__BANNER__", banner)
            .replace("__UNAVAILABLE__", unav).replace("__GSEATAB__", gsea_tab).replace("__GSEAPANEL__", gsea_panel))
    head_js = "window.__ttErrors=[];window.addEventListener('error',function(e){window.__ttErrors.push((e.message||'error')+(e.lineno?' @'+e.lineno:''));});"
    return ("<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">"
            f"<meta name=\"color-scheme\" content=\"light dark\"><title>{_esc(m['title'])}</title>\n<style>{CSS}</style>\n"
            f"<script>{head_js}</script>\n<script>{plotly_js}</script></head>\n<body>{body}\n"
            f"<script id=\"ttdata\" type=\"application/json\">{to_json(data)}</script>\n<script>{JS}</script>\n</body></html>\n")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def _parse_gsea(specs: Optional[List[str]]) -> Dict[str, Dict[str, str]]:
    out: Dict[str, Dict[str, str]] = {}
    for s in specs or []:
        parts = s.split("=", 2)
        if len(parts) != 3 or not all(parts):
            raise ValueError(f"--gsea expects COLLECTION=CONTRAST=PATH, got {s!r}")
        out.setdefault(parts[0], {})[parts[1]] = parts[2]
    return out


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Generate the static TTpipe results explorer.")
    ap.add_argument("--degs", nargs="+", required=True, help="DESeq2 result csv files (degs_<comparison>.csv)")
    ap.add_argument("--gene-info", default=None)
    ap.add_argument("--norm", default=None, help="normalized_counts.tsv")
    ap.add_argument("--samples", default=None, help="samples_deseq2.txt")
    ap.add_argument("--deg-summary", default=None)
    ap.add_argument("--gsea", action="append", default=None, help="COLLECTION=CONTRAST=PATH (repeatable)")
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--lfc", type=float, default=0.0)
    ap.add_argument("--condition-order", default="", help="comma-separated")
    ap.add_argument("--banner", default="")
    ap.add_argument("--title", default="TTpipe results explorer")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    degs = {contrast_name_from_path(p): p for p in a.degs}
    data = prepare_data(degs, gene_info=a.gene_info, norm=a.norm, samples=a.samples, gsea=_parse_gsea(a.gsea),
                        deg_summary=a.deg_summary, alpha=a.alpha, lfc=a.lfc,
                        condition_order=[c for c in a.condition_order.split(",") if c], banner=a.banner, title=a.title)
    doc = build_html(data)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(doc)
    mb = os.path.getsize(a.out) / 1e6
    print(f"results_explorer: wrote {a.out} ({mb:.1f} MB; {len(data['contrasts'])} contrast(s), "
          f"{data['meta']['n_genes']} genes, GSEA: {'yes' if data['gsea'] else 'no'}, "
          f"normalized counts: {'yes' if data['norm'] else 'no'})", file=sys.stderr)
    if mb > SIZE_WARN_MB:
        print(f"results_explorer: WARNING: the file is {mb:.1f} MB (> {SIZE_WARN_MB:.0f} MB); it may be slow to open. "
              "Fewer contrasts or GSEA collections make it smaller.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
