"""Tests for scripts/qc_summary.py and scripts/qc_report.py on SYNTHETIC data (tests/synth.py).

All data are invented; passing these tests shows the code parses the assumed file layouts and degrades
sensibly, not that the metrics are right for real runs.
"""
import os
import re
import sys

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))

import qc_report as q  # noqa: E402
import qc_summary as qs  # noqa: E402
import synth  # noqa: E402

CONFIG = os.path.join(HERE, "..", "config", "TTpipeConfig.yaml")


def run_all(run, tmp_path, name="r.html"):
    """qc_summary then qc_report on a synth run; returns (sample_qc DataFrame, html text)."""
    out = run["outdir"]
    tsv = os.path.join(out, "qc", "sample_qc.tsv")
    assert qs.main(["--outdir", out, "--samplesheet", run["samplesheet"], "--config", run["config"],
                    "--counts", run["counts"] or "", "--out", tsv]) == 0
    html_path = tmp_path / name
    args = ["--outdir", out, "--samplesheet", run["samplesheet"], "--config", run["config"], "--sample-qc", tsv,
            "--out", str(html_path), "--banner", synth.SYNTH_NOTICE]
    if run["counts"]:
        args += ["--counts", run["counts"]]
    assert q.main(args) == 0
    return pd.read_csv(tsv, sep="\t"), html_path.read_text()


def assert_self_contained(doc):
    # plotly.js itself contains URL strings; only markup outside <script> blocks can load external resources
    assert "Plotly.newPlot" in doc or "plotly" in doc.lower()
    inline_js = re.findall(r"<script>", doc)
    assert inline_js and not re.search(r"<script[^>]+\bsrc\s*=", doc, re.I)
    doc = re.sub(r"<script>.*?</script>", "", doc, flags=re.S)
    assert not re.search(r"<script[^>]+\bsrc\s*=", doc, re.I)
    assert not re.search(r"<link[^>]+\bhref\s*=", doc, re.I)
    assert not re.search(r"""(src|href)\s*=\s*["']\s*(https?:)?//""", doc, re.I)
    assert not re.search(r"@import|url\(\s*[\"']?https?:", doc)


# ---------------------------------------------------------------- parsers

def test_star_log_parser(tmp_path):
    p = tmp_path / "Log.final.out"
    p.write_text(synth.star_log_text(1_000_000, 80.0, 6.0, 1.0, 0.5, 8.0, 4.5))
    d = qs.parse_star_log(str(p))
    assert d["input_pairs"] == 1_000_000
    assert d["unique_pct"] == pytest.approx(80.0)
    assert d["multi_pct"] == pytest.approx(7.0)            # multiple + too many loci
    assert d["too_short_pct"] == pytest.approx(8.0)
    assert d["other_unmapped_pct"] == pytest.approx(5.0)   # mismatches + other


def test_star_log_garbage_gives_none(tmp_path):
    p = tmp_path / "x.out"
    p.write_text("nothing useful here\n")
    assert all(v is None for v in qs.parse_star_log(str(p)).values())


@pytest.mark.parametrize("comma", [False, True])
def test_picard_parser(tmp_path, comma):
    p = tmp_path / "m.txt"
    p.write_text(synth.picard_text(1000, 0.25, decimal_comma=comma))
    assert qs.parse_picard_metrics(str(p))["duplicate_pct"] == pytest.approx(25.0)


def test_picard_fallback_and_garbage(tmp_path):
    p = tmp_path / "m.txt"
    p.write_text("## METRICS CLASS\tx\nLIBRARY\tREAD_PAIRS_EXAMINED\tREAD_PAIR_DUPLICATES\nL\t200\t50\n")
    assert qs.parse_picard_metrics(str(p))["duplicate_pct"] == pytest.approx(25.0)   # derived, documented
    p.write_text("not a metrics file\n")
    assert qs.parse_picard_metrics(str(p))["duplicate_pct"] is None


def test_picard_two_libraries_weighted(tmp_path):
    p = tmp_path / "m.txt"
    p.write_text("## METRICS CLASS\tx\nLIBRARY\tREAD_PAIRS_EXAMINED\tPERCENT_DUPLICATION\nA\t100\t0.1\nB\t300\t0.5\n\n## HISTOGRAM\n")
    assert qs.parse_picard_metrics(str(p))["duplicate_pct"] == pytest.approx(40.0)


def test_pair_count(tmp_path):
    p = tmp_path / "c.txt"
    p.write_text("S_DMSO_rep1\t1234567\n")
    assert qs.parse_pair_count(str(p)) == 1234567
    p.write_text("")
    assert qs.parse_pair_count(str(p)) is None


@pytest.mark.parametrize("strandedness,expected", [("dUTP", 80.0), ("forward", 20.0)])
def test_star_genecounts_sense_depends_on_strandedness(tmp_path, strandedness, expected):
    p = tmp_path / "g.tab"
    p.write_text("N_unmapped\t9\t9\t9\nN_noFeature\t9\t9\t9\nG1\t100\t20\t80\nG2\t0\t0\t0\n")
    assert qs.parse_star_genecounts(str(p), strandedness)["sense_pct_star"] == pytest.approx(expected)


def test_star_genecounts_empty(tmp_path):
    p = tmp_path / "g.tab"
    p.write_text("N_unmapped\t9\t9\t9\n")
    assert qs.parse_star_genecounts(str(p), "dUTP")["sense_pct_star"] is None


# ---------------------------------------------------------------- count-matrix metrics

def _info():
    return pd.DataFrame({"gene_type": ["protein_coding", "rRNA", "protein_coding", "Mt_rRNA"],
                         "strand": ["+", "+", "-", "-"]}, index=pd.Index(["a", "b", "c", "d"], name="gene_id"))


def test_sense_pct_from_strand_matrices():
    fwd = pd.DataFrame({"s1": [90, 10, 5, 0]}, index=["a", "b", "c", "d"])
    rev = pd.DataFrame({"s1": [10, 0, 95, 10]}, index=["a", "b", "c", "d"])
    # sense: a fwd 90, b fwd 10, c rev 95, d rev 10 = 205 ; antisense = 10+0+5+0 = 15
    assert qs.sense_pct_from_strand_matrices(fwd, rev, _info())["s1"] == pytest.approx(100 * 205 / 220)


def test_rrna_pct_and_version_suffix_matching():
    m = pd.DataFrame({"s1": [70, 10, 10, 10]}, index=["a.1", "b.1", "c.2", "d.3"])
    assert qs.rrna_pct_from_matrix(m, _info())["s1"] == pytest.approx(20.0)   # b + d; ids matched without version
    assert qs.rrna_pct_from_matrix(pd.DataFrame({"s1": [1]}, index=["zzz"]), _info()) == {}


def test_replicate_correlation_values():
    rng = np.random.default_rng(0)
    base = rng.poisson(50, 300) + 1
    m = pd.DataFrame({"a1": base, "a2": base, "b1": base, "b2": rng.permutation(base)}, index=[f"g{i}" for i in range(300)])
    per, pairs = qs.replicate_correlations(qs.log_cpm(m), {"a1": "A", "a2": "A", "b1": "B", "b2": "B"})
    assert per["a1"] == pytest.approx(1.0) and per["a2"] == pytest.approx(1.0)
    assert per["b1"] < 0.5
    assert {p[0] for p in pairs} == {"A", "B"}
    # singleton condition -> no value
    per2, _ = qs.replicate_correlations(qs.log_cpm(m[["a1", "b1"]]), {"a1": "A", "b1": "B"})
    assert per2 == {}


# ---------------------------------------------------------------- traffic lights

@pytest.mark.parametrize("value,direction,expected", [
    (85, "higher_is_better", "pass"), (80, "higher_is_better", "pass"), (79.9, "higher_is_better", "warn"),
    (60, "higher_is_better", "warn"), (59.9, "higher_is_better", "fail"),
    (30, "lower_is_better", "pass"), (40, "lower_is_better", "pass"), (41, "lower_is_better", "warn"),
    (60, "lower_is_better", "warn"), (61, "lower_is_better", "fail"), (None, "lower_is_better", "na"),
    (float("nan"), "higher_is_better", "na"),
])
def test_classify(value, direction, expected):
    warn, fail = (80, 60) if direction == "higher_is_better" else (40, 60)
    assert q.classify(value, {"direction": direction, "warn": warn, "fail": fail}) == expected


def test_overall_is_worst_and_na_ignored():
    thr = {"a": {"direction": "higher_is_better", "warn": 10, "fail": 5},
           "b": {"direction": "lower_is_better", "warn": 10, "fail": 20}}
    assert q.classify_sample({"a": 20, "b": 1}, thr)["overall"] == "pass"
    assert q.classify_sample({"a": 7, "b": 1}, thr)["overall"] == "warn"
    assert q.classify_sample({"a": 7, "b": 50}, thr)["overall"] == "fail"
    assert q.classify_sample({"a": None, "b": 1}, thr)["overall"] == "pass"
    assert q.classify_sample({}, thr)["overall"] == "na"


def test_thresholds_from_repo_config_match_defaults():
    thr = q.load_thresholds(qs.load_config(CONFIG))
    assert set(thr) == set(q.DEFAULT_THRESHOLDS)
    assert thr["sense_pct"] == {"direction": "higher_is_better", "warn": 90, "fail": 75}
    assert q.load_thresholds({}) == q.DEFAULT_THRESHOLDS


def test_pca_shapes_and_separation():
    rng = np.random.default_rng(1)
    base = rng.normal(5, 1, 200)
    cols = {}
    for i, grp in enumerate(["x", "x", "y", "y"]):
        v = base + rng.normal(0, 0.05, 200)
        if grp == "y":
            v[:40] += 3
        cols[f"s{i}"] = v
    coords, expl = q.compute_pca(pd.DataFrame(cols))
    assert coords.shape[0] == 4 and expl[0] > 0.8
    assert np.sign(coords.loc["s0", "PC1"]) == np.sign(coords.loc["s1", "PC1"]) != np.sign(coords.loc["s2", "PC1"])
    assert q.compute_pca(pd.DataFrame(cols).iloc[:, :2]) is None


# ---------------------------------------------------------------- end to end (synthetic)

SECTIONS = ["run", "skipped", "summary", "alignment", "dup", "strand", "rrna", "correlation", "pca", "deg",
            "tracks", "notassessed", "methods"]


def test_full_run_report(tmp_path):
    run = synth.make_run(str(tmp_path / "out"), "full", poor=("S_2H_rep2",))
    qc, doc = run_all(run, tmp_path)
    by = qc.set_index("sample")
    assert by.loc["S_DMSO_rep1", "unique_pct"] == pytest.approx(88.0)
    assert by.loc["S_DMSO_rep1", "duplicate_pct"] == pytest.approx(25.0)
    assert by.loc["S_DMSO_rep1", "input_pairs_million"] == pytest.approx(40.0)
    assert by.loc["S_DMSO_rep1", "sense_pct_htseq"] == pytest.approx(97.0, abs=1.0)
    assert by.loc["S_DMSO_rep1", "sense_pct_star"] == pytest.approx(97.0, abs=1.0)
    assert by.loc["S_DMSO_rep1", "sense_source"] == "htseq"
    assert by.loc["S_DMSO_rep1", "rrna_pct"] > 0
    assert by.loc["S_DMSO_rep1", "replicate_correlation"] > 0.9
    assert by.loc["S_2H_rep2", "sense_pct"] == pytest.approx(55.0, abs=2.0)
    # report
    for sec in SECTIONS:
        assert f'id="{sec}"' in doc, sec
    assert synth.SYNTH_NOTICE in doc
    assert_self_contained(doc)
    for pid in ("plot-align", "plot-corr", "plot-pca", "plot-deg", "plot-strand", "plot-dup", "plot-rrna", "plot-repcorr"):
        assert pid in doc, pid
    assert "FAIL" in doc and "PASS" in doc
    skipped = doc.split('id="skipped"')[1].split("</section>")[0]
    assert "None: every section had data" in skipped
    assert "VST counts from DESeq2" in doc
    na = doc.split('id="notassessed"')[1].split("</section>")[0]
    for term in ("Spike-in fraction", "Labeled-vs-total enrichment", "Synthesis and degradation"):
        assert term in na
    assert "S_DMSO_rep1.bw" not in doc and "DMSO_forward.bw" in doc and "scale_factor" in doc
    assert "S_2H_rep2" in doc


def test_poor_sample_status(tmp_path):
    run = synth.make_run(str(tmp_path / "out"), "full", poor=("S_2H_rep2",))
    qc, _ = run_all(run, tmp_path)
    thr = q.load_thresholds(qs.load_config(run["config"]))
    st = {r["sample"]: q.classify_sample({k: r[k] for k in thr}, thr) for _, r in qc.iterrows()}
    assert st["S_2H_rep2"]["overall"] == "fail"
    assert st["S_2H_rep2"]["unique_pct"] == "fail" and st["S_2H_rep2"]["sense_pct"] == "fail"
    assert st["S_DMSO_rep1"]["overall"] in ("pass", "warn")


def test_processing_only_degradation(tmp_path):
    run = synth.make_run(str(tmp_path / "out"), "processing")
    _qc, doc = run_all(run, tmp_path)
    skipped = doc.split('id="skipped"')[1].split("</section>")[0]
    assert "Differential expression" in skipped
    assert "computed here from" in doc          # PCA falls back to log-CPM
    assert "plot-pca" in doc and "plot-deg" not in doc
    assert "plot-align" in doc
    assert_self_contained(doc)


def test_analysis_only_degradation(tmp_path):
    run = synth.make_run(str(tmp_path / "out"), "analysis")
    # analysis-only: no per-sample stats exist and no processing flags are given by the config
    out = run["outdir"]
    tsv = os.path.join(out, "qc", "sample_qc.tsv")
    assert qs.main(["--outdir", out, "--samplesheet", run["samplesheet"], "--config", run["config"],
                    "--counts", run["counts"], "--out", tsv]) == 0
    qc = pd.read_csv(tsv, sep="\t")
    for col in ("unique_pct", "duplicate_pct", "sense_pct", "input_pairs_million"):
        assert qc[col].isna().all(), col       # blank, never fabricated
    assert qc["replicate_correlation"].notna().all()
    assert qc["rrna_pct"].notna().all()
    html_path = tmp_path / "a.html"
    assert q.main(["--outdir", out, "--samplesheet", run["samplesheet"], "--config", run["config"], "--sample-qc", tsv,
                   "--counts", run["counts"], "--run-processing", "0", "--run-analysis", "1", "--out", str(html_path),
                   "--banner", synth.SYNTH_NOTICE]) == 0
    doc = html_path.read_text()
    skipped = doc.split('id="skipped"')[1].split("</section>")[0]
    for what in ("Alignment", "Duplication", "Strandedness", "Tracks", "Differential expression"):
        assert what in skipped, what
    assert "run_processing is false" in skipped
    assert "plot-corr" in doc and "plot-pca" in doc and "plot-rrna" in doc and "plot-repcorr" in doc
    assert "plot-align" not in doc and "plot-strand" not in doc and "plot-dup" not in doc
    assert synth.SYNTH_NOTICE in doc
    assert_self_contained(doc)


def test_report_without_sample_qc_file_recomputes(tmp_path):
    run = synth.make_run(str(tmp_path / "out"), "full")
    html_path = tmp_path / "r.html"
    assert q.main(["--outdir", run["outdir"], "--samplesheet", run["samplesheet"], "--config", run["config"],
                   "--out", str(html_path)]) == 0
    assert "plot-align" in html_path.read_text()


def test_decimal_comma_picard_end_to_end(tmp_path):
    run = synth.make_run(str(tmp_path / "out"), "full", decimal_comma=True)
    qc, _ = run_all(run, tmp_path)
    assert qc["duplicate_pct"].iloc[0] == pytest.approx(25.0)


def test_empty_outdir_still_writes_a_report(tmp_path):
    sheet = tmp_path / "s.txt"
    sheet.write_text("Name\tCondition\tReplicate\nA\tc1\tr1\nB\tc1\tr2\n")
    html_path = tmp_path / "x" / "r.html"
    assert q.main(["--outdir", str(tmp_path / "empty"), "--samplesheet", str(sheet), "--out", str(html_path)]) == 0
    doc = html_path.read_text()
    assert "Skipped" in doc and "No QC data" in doc
    assert_self_contained(doc)


def test_sample_names_are_escaped(tmp_path):
    run = synth.make_run(str(tmp_path / "out"), "processing")
    # only the report-side escaping is under test: put the odd name in the qc table
    tsv = os.path.join(run["outdir"], "qc", "sample_qc.tsv")
    qs.main(["--outdir", run["outdir"], "--samplesheet", run["samplesheet"], "--config", run["config"], "--out", tsv])
    t = pd.read_csv(tsv, sep="\t")
    t.loc[0, "sample"] = "a<b>&c"
    t.to_csv(tsv, sep="\t", index=False)
    out = tmp_path / "e.html"
    q.main(["--outdir", run["outdir"], "--samplesheet", run["samplesheet"], "--config", run["config"],
            "--sample-qc", tsv, "--out", str(out)])
    doc = out.read_text()
    assert "a&lt;b&gt;&amp;c" in doc and "<b>&c" not in doc
