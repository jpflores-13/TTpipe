"""Tests for scripts/results_explorer.py on SYNTHETIC tables built in this file.

Every table below is INVENTED by a random number generator (gene IDs SYN00000..., contrasts of fake
conditions). Passing these tests shows that the data preparation, the HTML and the JavaScript behave on the
assumed file layouts (DESeq2 csv, gene_info.tsv, normalized_counts.tsv, fgsea tsv); it says nothing about
real TT-seq results. The browser test runs the page in headless Chrome when Chrome is installed.
"""
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))

import results_explorer as rx  # noqa: E402

N = 300


def make_tables(root, contrasts=("2H_vs_DMSO", "6H_vs_DMSO"), norm=True, gsea=False, gene_info=True, seed=3):
    """Write synthetic DESeq2 / gene / count / GSEA tables under `root`; return a dict of paths."""
    rng = np.random.default_rng(seed)
    root = Path(root)
    d = root / "analysis" / "deseq2"
    d.mkdir(parents=True, exist_ok=True)
    ids = [f"SYN{i:05d}" for i in range(N)]
    paths = {"degs": {}, "gene_info": None, "norm": None, "samples": None, "summary": None, "gsea": {}}
    for ci, name in enumerate(contrasts):
        base = rng.lognormal(5, 1.5, N)
        lfc = rng.normal(0, 0.6, N)
        lfc[:40] += 2.5 * (1 if ci == 0 else -1)
        padj = np.clip(10 ** (-np.abs(lfc) * 4 + rng.normal(0, 0.5, N)), 1e-300, 1)
        df = pd.DataFrame({"gene_id": ids, "baseMean": base, "log2FoldChange": lfc, "lfcSE": 0.2, "stat": lfc / 0.2,
                           "pvalue": padj, "padj": padj})
        df.loc[5:9, "padj"] = np.nan            # NA padj, finite baseMean (independent filtering)
        df.loc[10:11, ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]] = np.nan   # all-NA rows
        p = d / f"degs_{name}.csv"
        df.to_csv(p, index=False)
        paths["degs"][name] = str(p)
    (root / "counts").mkdir(exist_ok=True)
    if gene_info:
        types = ["protein_coding"] * 150 + ["lncRNA"] * 80 + ["miRNA"] * 30 + ["snRNA"] * 20 + ["pseudogene"] * 20
        names = [f"Gene{i}" for i in range(N)]
        names[20] = ""                          # missing symbol -> falls back to gene_id
        pd.DataFrame({"gene_id": ids, "gene_name": names, "gene_type": types, "chr": "chrSyn", "start": 1, "end": 100,
                      "strand": "+", "length": 100}).to_csv(root / "counts" / "gene_info.tsv", sep="\t", index=False)
        paths["gene_info"] = str(root / "counts" / "gene_info.tsv")
    samples = [("S_DMSO_1", "DMSO"), ("S_DMSO_2", "DMSO"), ("S_2H_1", "2H"), ("S_2H_2", "2H"), ("S_6H_1", "6H"), ("S_6H_2", "6H")]
    pd.DataFrame({"sample": [s for s, _ in samples], "condition": [c for _, c in samples], "replicate": "r"}).to_csv(
        root / "counts" / "samples_deseq2.txt", sep="\t", index=False)
    paths["samples"] = str(root / "counts" / "samples_deseq2.txt")
    if norm:
        m = pd.DataFrame(rng.lognormal(5, 1, (N, len(samples))), index=ids, columns=[s for s, _ in samples])
        m.index.name = "gene_id"
        m.to_csv(d / "normalized_counts.tsv", sep="\t")
        paths["norm"] = str(d / "normalized_counts.tsv")
    pd.DataFrame({"comparison": list(contrasts), "n_up": 30, "n_down": 5, "n_tested": N - 2}).to_csv(
        d / "deg_summary.tsv", sep="\t", index=False)
    paths["summary"] = str(d / "deg_summary.tsv")
    if gsea:
        g = root / "analysis" / "gsea"
        g.mkdir(parents=True, exist_ok=True)
        for name in contrasts:
            rows = [(f"HALLMARK_SYN_{i}", 10 + i, 0.5, (-1) ** i * 1.5, 10 ** -(i + 1), None if i == 7 else 10 ** -(i + 1) * 2,
                     ";".join(f"Gene{j}" for j in range(i * 5, i * 5 + 30))) for i in range(10)]
            p = g / f"Hallmarks_{name}.tsv"
            pd.DataFrame(rows, columns=["pathway", "size", "ES", "NES", "pval", "padj", "leadingEdge"]).to_csv(p, sep="\t", index=False)
            paths["gsea"].setdefault("Hallmarks", {})[name] = str(p)
    return paths


def build(paths, **kw):
    data = rx.prepare_data(paths["degs"], gene_info=paths["gene_info"], norm=paths["norm"], samples=paths["samples"],
                           gsea=paths["gsea"] or None, deg_summary=paths["summary"], alpha=0.01, lfc=0.0,
                           banner="SYNTHETIC TEST DATA", **kw)
    return data, rx.build_html(data)


def embedded_json(doc):
    m = re.search(r'<script id="ttdata" type="application/json">(.*?)</script>', doc, re.S)
    assert m
    return json.loads(m.group(1))


# ---- data preparation ----

def test_sig_round_and_json_safety():
    assert rx.sig_round(1.23456789, 5) == 1.2346
    assert rx.sig_round(0.000123456789, 3) == 0.000123
    assert rx.sig_round(float("nan")) is None and rx.sig_round(None) is None and rx.sig_round(float("inf")) is None
    assert rx.sig_round(0) == 0.0
    s = rx.to_json({"a": "</script><b>&"})
    assert "</script>" not in s and "<" not in s and json.loads(s)["a"] == "</script><b>&"
    with pytest.raises(ValueError):
        rx.to_json({"a": float("nan")})


def test_contrast_name_from_path():
    assert rx.contrast_name_from_path("/x/degs_2H_vs_DMSO.csv") == "2H_vs_DMSO"
    assert rx.contrast_name_from_path("/x/other.csv") == "other"


def test_prepare_multi_contrast(tmp_path):
    p = make_tables(tmp_path)
    data, _ = build(p)
    assert [c["name"] for c in data["contrasts"]] == ["2H_vs_DMSO", "6H_vs_DMSO"]
    assert data["meta"]["n_genes"] == N - 2                      # the two all-NA genes are dropped
    for c in data["contrasts"]:
        assert len(c["g"]) == len(c["bm"]) == len(c["lfc"]) == len(c["padj"]) == N - 2
    g = data["genes"]
    assert len(g["id"]) == len(g["name"]) == len(g["type"]) == N - 2
    assert "SYN00010" not in g["id"] and "SYN00011" not in g["id"]
    assert data["meta"]["conditions"] == ["DMSO", "2H", "6H"]
    assert data["meta"]["samples"][0] == "S_DMSO_1"
    assert data["meta"]["deg_summary"][0][0] == "2H_vs_DMSO"
    assert data["gsea"] == {}


def test_na_padj_kept_as_null_and_values_rounded(tmp_path):
    p = make_tables(tmp_path)
    data, doc = build(p)
    c = data["contrasts"][0]
    pos = {g: k for k, g in enumerate(data["genes"]["id"][i] for i in c["g"])}
    for gid in ("SYN00005", "SYN00009"):
        k = pos[gid]
        assert c["padj"][k] is None and c["bm"][k] is not None       # NA padj stays, shown as n.s. by the page
    assert "NaN" not in doc.split('id="ttdata"')[1].split("</script>")[0]
    for v in c["bm"]:
        assert len(f"{v:.10g}".replace(".", "").replace("-", "").lstrip("0").split("e")[0]) <= 6 or v == 0


def test_biotype_order_and_gene_name_fallback(tmp_path):
    p = make_tables(tmp_path)
    data, _ = build(p)
    assert data["biotypes"][:2] == ["protein_coding", "lncRNA"]
    assert set(data["biotypes"]) == {"protein_coding", "lncRNA", "miRNA", "snRNA", "pseudogene"}
    assert data["biotypes"][2] == "miRNA"                         # others by decreasing count (miRNA 30 > snRNA 20 = pseudogene 20)
    gi = data["genes"]["id"].index("SYN00020")
    assert data["genes"]["name"][gi] == "SYN00020"                # empty gene_name -> gene_id
    assert data["meta"]["n_unannotated"] == 0


def test_missing_gene_info_degrades(tmp_path):
    p = make_tables(tmp_path, gene_info=False)
    data, doc = build(p)
    assert data["biotypes"] == ["unknown"]
    assert data["genes"]["name"] == data["genes"]["id"]
    assert any("gene_info" in u for u in data["meta"]["unavailable"])
    assert "Not available in this file" in doc


def test_missing_norm_counts(tmp_path):
    p = make_tables(tmp_path, norm=False)
    data, doc = build(p)
    assert data["norm"] is None
    assert any("normalized_counts" in u for u in data["meta"]["unavailable"])
    assert "Normalized counts (normalized_counts.tsv) were not found" in doc


def test_norm_block_aligned_with_genes(tmp_path):
    p = make_tables(tmp_path)
    data, _ = build(p)
    nb = data["norm"]
    assert nb["samples"] == data["meta"]["samples"] and nb["cond"][0] == "DMSO" and nb["cond"][-1] == "6H"
    assert len(nb["v"]) == data["meta"]["n_genes"] and len(nb["v"][0]) == 6
    raw = pd.read_csv(p["norm"], sep="\t", index_col=0)
    gid = data["genes"]["id"][7]
    assert nb["v"][7][2] == pytest.approx(raw.loc[gid, "S_2H_1"], rel=1e-3)


def test_gsea_sorted_and_tab(tmp_path):
    p = make_tables(tmp_path, gsea=True)
    data, doc = build(p)
    rows = data["gsea"]["Hallmarks"]["2H_vs_DMSO"]
    padj = [r[4] for r in rows]
    assert padj[-1] is None and padj[:-1] == sorted(x for x in padj if x is not None)    # NA last, ascending
    assert data["gsea"]["Hallmarks"].keys() == {"2H_vs_DMSO", "6H_vs_DMSO"}
    assert 'data-panel="p-gsea"' in doc and 'id="p-gsea"' in doc


def test_no_gsea_means_no_tab(tmp_path):
    p = make_tables(tmp_path)
    _, doc = build(p)
    body = re.sub(r"<script.*?</script>", "", doc, flags=re.S)
    assert "p-gsea" not in body and ">GSEA<" not in body
    assert "no GSEA tab" in body
    # a collection whose files are all missing also gives no tab
    gs = {"H": {"2H_vs_DMSO": str(tmp_path / "nope.tsv")}}
    data = rx.prepare_data(p["degs"], gsea=gs)
    assert data["gsea"] == {}


def test_single_contrast(tmp_path):
    p = make_tables(tmp_path, contrasts=("2H_vs_DMSO",))
    data, doc = build(p)
    assert len(data["contrasts"]) == 1
    assert embedded_json(doc)["contrasts"][0]["name"] == "2H_vs_DMSO"


def test_no_deg_files_is_an_error():
    with pytest.raises(ValueError):
        rx.prepare_data({})


# ---- HTML ----

def strip_scripts(doc):
    return re.sub(r"<script\b.*?</script>", "", doc, flags=re.S)


def test_html_self_contained_and_json_parses(tmp_path):
    p = make_tables(tmp_path, gsea=True)
    data, doc = build(p)
    assert len(doc) > 1_000_000 and "plotly.js" in doc and doc.lower().count("<script>") + doc.lower().count("<script ") >= 3
    assert not re.search(r"<script[^>]+\bsrc\s*=", doc, re.I)
    body = strip_scripts(doc)
    assert not re.search(r"<link[^>]+\bhref\s*=", body, re.I)
    assert not re.search(r"""(src|href)\s*=\s*["']\s*(https?:)?//""", body, re.I)
    assert not re.search(r"@import|url\(\s*[\"']?https?:", doc.split("</style>")[0])
    # no network calls in OUR javascript (the plotly bundle is excluded)
    ours = doc.split('<script id="ttdata"')[1]
    assert not re.search(r"fetch\(|XMLHttpRequest|WebSocket|https?://", ours.split("<script>")[1])
    emb = embedded_json(doc)
    assert [c["name"] for c in emb["contrasts"]] == ["2H_vs_DMSO", "6H_vs_DMSO"]
    assert "SYN00000" in emb["genes"]["id"] and "SYN00010" not in emb["genes"]["id"]
    assert emb["meta"]["banner"] == "SYNTHETIC TEST DATA" and "SYNTHETIC TEST DATA" in body
    assert emb["meta"]["alpha"] == 0.01


def test_cli_writes_file_and_logs_size(tmp_path, capsys):
    p = make_tables(tmp_path, gsea=True)
    out = tmp_path / "analysis" / "TTpipe_results_explorer.html"
    args = ["--degs", *p["degs"].values(), "--gene-info", p["gene_info"], "--norm", p["norm"], "--samples", p["samples"],
            "--deg-summary", p["summary"], "--alpha", "0.05", "--lfc", "1", "--condition-order", "DMSO,2H,6H",
            "--out", str(out)]
    for coll, per in p["gsea"].items():
        for cname, path in per.items():
            args += ["--gsea", f"{coll}={cname}={path}"]
    assert rx.main(args) == 0
    err = capsys.readouterr().err
    assert "MB" in err and "WARNING" not in err
    emb = embedded_json(out.read_text())
    assert emb["meta"]["alpha"] == 0.05 and emb["meta"]["lfc"] == 1.0 and emb["gsea"]["Hallmarks"]
    with pytest.raises(ValueError):
        rx.main(["--degs", p["degs"]["2H_vs_DMSO"], "--gsea", "bad", "--out", str(out)])


def test_size_warning_threshold(tmp_path, capsys, monkeypatch):
    p = make_tables(tmp_path, contrasts=("a_vs_b",))
    monkeypatch.setattr(rx, "SIZE_WARN_MB", 0.0001)
    assert rx.main(["--degs", *p["degs"].values(), "--out", str(tmp_path / "x.html")]) == 0
    assert "WARNING" in capsys.readouterr().err


def test_argparse_accepts_degs_and_no_choices():
    # regression guard for the Python 3.11 argparse trap (choices= together with nargs="*" defaults)
    src = (Path(HERE).parent / "scripts" / "results_explorer.py").read_text()
    assert "choices=" not in src


# ---- the JavaScript really runs (headless Chrome, skipped when Chrome is absent) ----

def find_chrome():
    cands = ["/Applications/Google Chrome.app/Contents/MacOS/Google Chrome", "/Applications/Chromium.app/Contents/MacOS/Chromium"]
    for n in ("google-chrome", "google-chrome-stable", "chromium", "chromium-browser", "chrome"):
        w = shutil.which(n)
        if w:
            cands.insert(0, w)
    for c in cands:
        if os.path.exists(c):
            return c
    return None


CHROME = find_chrome()
needs_chrome = pytest.mark.skipif(CHROME is None, reason="Chrome/Chromium not installed")


def dump_dom(html_path, hash_=""):
    prof = Path(html_path).parent / "chrome_profile"
    cmd = [CHROME, "--headless=new", "--disable-gpu", "--no-sandbox", "--no-first-run", f"--user-data-dir={prof}",
           "--allow-file-access-from-files", "--virtual-time-budget=20000", "--dump-dom", f"file://{html_path}{hash_}"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        pytest.skip("headless Chrome timed out")
    if r.returncode != 0 or "<html" not in r.stdout:
        pytest.skip(f"headless Chrome unusable here (exit {r.returncode}): {r.stderr[-200:]}")
    return r.stdout


@needs_chrome
@pytest.mark.parametrize("variant", ["full", "no_norm_no_gsea", "single"])
def test_javascript_runs_in_headless_chrome(tmp_path, variant):
    kw = {"full": dict(gsea=True), "no_norm_no_gsea": dict(norm=False), "single": dict(contrasts=("2H_vs_DMSO",), gsea=True)}[variant]
    p = make_tables(tmp_path, **kw)
    _, doc = build(p)
    f = tmp_path / "explorer.html"
    f.write_text(doc, encoding="utf-8")
    # 1) plain load: controls rendered, plot and table drawn, no script error
    dom = dump_dom(f)
    assert 'data-ready="1"' in dom, "the page script did not finish"
    assert 'id="errbox" class="err" role="alert" hidden' in dom
    assert len(re.findall(r"<option value=\"\d+\">", dom)) == len(embedded_json(doc)["contrasts"])
    assert 'class="js-plotly-plot"' in dom or "js-plotly-plot" in dom
    assert "Up: " in dom and "genes shown" in dom
    assert len(re.findall(r'<tr data-gi="', dom)) > 0
    # 2) self-test through the real controls (#selftest hash): every check must PASS
    dom = dump_dom(f, "#selftest")
    m = re.search(r'<pre id="selftest">(.*?)</pre>', dom, re.S)
    assert m, "self-test output missing"
    lines = m.group(1).strip().splitlines()
    assert lines[-1] == "SELFTEST DONE failures=0", "\n".join(lines)
    assert not [ln for ln in lines if ln.startswith("FAIL")], "\n".join(lines)
    assert sum(ln.startswith("PASS") for ln in lines) >= 15
    if variant == "no_norm_no_gsea":
        assert any("count plot unavailable message" in ln for ln in lines)
        assert any("no GSEA tab" in ln for ln in lines)
    if variant == "full":
        assert any("GSEA tab visible" in ln for ln in lines) and any("second contrast" in ln for ln in lines)
