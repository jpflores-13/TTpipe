"""Snakemake dry runs of TTpipe against the PLACEHOLDER fixtures in tests/dryrun/ (see its README).

A dry run (-n) parses the samplesheet and config, builds the DAG and checks that input files exist; it
runs no tools. Skipped when snakemake is not installed. Run from anywhere: commands run from the repo root.
Fixture: 2 conditions (DMSO, 2H) x 2 replicates = 4 samples, so e.g. 4 star_align jobs.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MAIN = "config/TTpipeConfig.yaml"
FULL = "tests/dryrun/config.yaml"
PROC = "tests/dryrun/config_processing_only.yaml"
ANA = "tests/dryrun/config_analysis_only.yaml"

# also find snakemake next to the running interpreter (a venv used without being activated)
SNAKEMAKE = shutil.which("snakemake", path=os.pathsep.join([os.path.dirname(sys.executable), os.environ.get("PATH", "")]))
needs_snakemake = pytest.mark.skipif(SNAKEMAKE is None, reason="snakemake not installed")


def dry_run(*configs, **kv):
    """Dry run with the main config overlaid by `configs`; kv become --config key=value (use 'True'/'False')."""
    cmd = [SNAKEMAKE, "-n", "-s", "workflows/TTpipe.snakefile", "--cores", "1", "--configfile", MAIN, *configs]
    if kv:
        cmd += ["--config", *[f"{k}={v}" for k, v in kv.items()]]
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)


def jobs(r):
    """Rule -> job count from the first 'Job stats' table of a dry run."""
    out, started = {}, False
    for ln in r.stdout.splitlines():
        if ln.startswith("Job stats"):
            started = True
            continue
        if started:
            p = ln.split()
            if len(p) == 2 and p[1].isdigit():
                out[p[0]] = int(p[1])
            elif not ln.strip() and out:
                break
    return out


@needs_snakemake
def test_full_run():
    r = dry_run(FULL)
    assert r.returncode == 0, r.stderr
    j = jobs(r)
    assert j["star_align"] == 4 and j["fastqc"] == 8 and j["deseq2"] == 1
    assert j["qc_summary"] == 1 and j["qc_report"] == 1
    assert j["bigwig"] == 6            # 2 conditions x (all, forward, reverse)
    assert "fastp" not in j and "bigwig_sample" not in j and "fgsea" not in j
    assert "output_ci/qc/TTpipe_qc_report.html" in r.stdout


@needs_snakemake
def test_processing_only():
    r = dry_run(PROC)
    assert r.returncode == 0, r.stderr
    j = jobs(r)
    assert j["star_align"] == 4 and j["qc_report"] == 1
    assert "deseq2" not in j and "sankey" not in j and "sample_table" not in j


@needs_snakemake
def test_analysis_only_has_no_alignment():
    r = dry_run(ANA)
    assert r.returncode == 0, r.stderr
    j = jobs(r)
    assert "star_align" not in j and "fastqc" not in j and "htseq_count" not in j
    assert j["deseq2"] == 1 and j["qc_summary"] == 1 and j["qc_report"] == 1
    assert j["results_explorer"] == 1
    assert set(j) - {"multiqc"} == {"all", "deseq2", "gene_info", "qc_summary", "qc_report", "results_explorer",
                                    "sample_table", "sankey", "total"}


@needs_snakemake
def test_trim_adds_fastp():
    base, j = jobs(dry_run(FULL)), jobs(dry_run(FULL, trim="True"))
    assert j["fastp"] == 4 and j["total"] == base["total"] + 4


@needs_snakemake
def test_sample_tracks_add_per_sample_bigwigs():
    base, j = jobs(dry_run(FULL)), jobs(dry_run(FULL, sample_tracks="True"))
    assert j["bigwig_sample"] == 12    # 4 samples x (all, forward, reverse)
    assert j["total"] == base["total"] + 12


@needs_snakemake
def test_gsea_collections_add_fgsea():
    r = dry_run(FULL, gsea_collections="{H: tests/dryrun/placeholder.gmt}")
    assert r.returncode == 0, r.stderr
    j = jobs(r)
    assert j["fgsea"] == 1 and j["rnk"] == 1       # 1 collection x 1 comparison (2H_vs_DMSO)
    assert "gsea/H_2H_vs_DMSO.tsv" in r.stdout


@needs_snakemake
def test_qc_report_switch_removes_report_jobs():
    base, j = jobs(dry_run(FULL)), jobs(dry_run(FULL, qc_report="False"))
    assert "qc_report" not in j and "qc_summary" not in j and j["total"] == base["total"] - 2


@needs_snakemake
def test_results_explorer_switch_adds_or_removes_one_job():
    on, off = jobs(dry_run(FULL, results_explorer="True")), jobs(dry_run(FULL, results_explorer="False"))
    assert on["results_explorer"] == 1 and "results_explorer" not in off
    assert on["total"] == off["total"] + 1
    assert jobs(dry_run(FULL)) == on               # shipped config: on by default
    r = dry_run(FULL)
    assert "output_ci/analysis/TTpipe_results_explorer.html" in r.stdout


@needs_snakemake
def test_results_explorer_in_analysis_only_mode_and_with_gsea():
    on, off = jobs(dry_run(ANA, results_explorer="True")), jobs(dry_run(ANA, results_explorer="False"))
    assert on["results_explorer"] == 1 and "results_explorer" not in off and on["total"] == off["total"] + 1
    r = dry_run(ANA, gsea_collections="{H: tests/dryrun/placeholder.gmt}")
    assert r.returncode == 0, r.stderr
    j = jobs(r)
    assert j["results_explorer"] == 1 and j["fgsea"] == 1         # still exactly one explorer job, after fgsea
    block = r.stdout.split("rule results_explorer:")[1].split("\n\n")[0]
    assert "gsea/H_2H_vs_DMSO.tsv" in block and "degs_2H_vs_DMSO.csv" in block and "normalized_counts.tsv" in block


@needs_snakemake
def test_results_explorer_not_in_processing_only_mode():
    assert "results_explorer" not in jobs(dry_run(PROC))


@needs_snakemake
@pytest.mark.parametrize("kv,msg", [
    ({"strandedness": "bad"}, "'strandedness' must be 'dUTP' or 'forward'"),
    ({"count_strand": "bad"}, "'count_strand' must be 'sum' or 'sense'"),
    ({"reference_condition": "nope"}, "reference_condition 'nope' is not a Condition"),
])
def test_bad_config_is_rejected(kv, msg):
    r = dry_run(FULL, **kv)
    assert r.returncode != 0
    assert msg in (r.stderr + r.stdout)


@needs_snakemake
def test_bad_samplesheets_are_rejected(tmp_path):
    good = (ROOT / "tests/dryrun/samplesheet.txt").read_text()
    cases = {
        "space": (good.replace("CI_DMSO_rep1", "bad name"), "may only contain letters, digits"),
        "dup": (good.replace("CI_DMSO_rep2", "CI_DMSO_rep1"), "duplicate Name"),
        "cols": ("\n".join("\t".join(ln.split("\t")[:2]) for ln in good.splitlines()) + "\n", "is missing column(s)"),
    }
    for name, (text, msg) in cases.items():
        sheet = tmp_path / f"{name}.txt"
        sheet.write_text(text)
        r = dry_run(FULL, samplesheet=str(sheet))
        assert r.returncode != 0, name
        assert msg in (r.stderr + r.stdout), name


@pytest.mark.skipif(shutil.which("bash") is None, reason="bash not installed")
@pytest.mark.parametrize("script", ["TTpipe.sh", "unlock.sh"])
def test_launcher_scripts_parse(script):
    r = subprocess.run(["bash", "-n", script], cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


# ---- stranded signal page (rule stranded_signal, workflows/rules/stranded.smk) ----

@needs_snakemake
def test_stranded_plots_switch_adds_or_removes_exactly_one_job():
    on, off = jobs(dry_run(FULL, stranded_plots="True")), jobs(dry_run(FULL, stranded_plots="False"))
    assert on["stranded_signal"] == 1 and "stranded_signal" not in off
    assert on["total"] == off["total"] + 1
    assert jobs(dry_run(FULL)) == on                       # shipped config: on by default
    assert "output_ci/analysis/TTpipe_stranded_signal.html" in dry_run(FULL).stdout


@needs_snakemake
def test_stranded_signal_inputs_follow_the_run_mode():
    r = dry_run(FULL)
    block = r.stdout.split("rule stranded_signal:")[1].split("\n\n")[0]
    for part in ("DMSO_forward.bw", "2H_reverse.bw", "scale_factors.tsv", "gene_info.tsv", "genes.gtf", "counts_forward.txt",
                 "counts_reverse.txt", "counts_deseq2.txt", "samplesheet.txt", "degs_2H_vs_DMSO.csv"):
        assert part in block, part
    r = dry_run(FULL, track_top_degs="0")
    assert "degs_2H_vs_DMSO.csv" not in r.stdout.split("rule stranded_signal:")[1].split("\n\n")[0]
    r = dry_run(PROC)                                      # processing only: no DESeq2 tables exist, none are requested
    assert r.returncode == 0, r.stderr
    j = jobs(r)
    assert j["stranded_signal"] == 1 and "deseq2" not in j
    assert "degs_" not in r.stdout.split("rule stranded_signal:")[1].split("\n\n")[0]


@needs_snakemake
def test_stranded_signal_skipped_cleanly_in_analysis_only_mode():
    r = dry_run(ANA)
    assert r.returncode == 0, r.stderr                     # no bigWigs: skipped, not an error
    assert "stranded_signal" not in jobs(r)
    r = dry_run(ANA, stranded_plots="True")
    assert r.returncode == 0 and "stranded_signal" not in jobs(r)


@needs_snakemake
def test_stranded_signal_in_analysis_only_when_bigwigs_already_exist(tmp_path):
    out = tmp_path / "out"
    for sub, names in {"tracks": ["DMSO_forward.bw", "DMSO_reverse.bw", "2H_forward.bw", "2H_reverse.bw", "scale_factors.tsv"],
                       "counts": ["gene_info.tsv", "counts_forward.txt", "counts_reverse.txt"]}.items():
        (out / sub).mkdir(parents=True)
        for n in names:
            (out / sub / n).write_text("")
    r = dry_run(ANA, outputDir=str(out))
    assert r.returncode == 0, r.stderr
    j = jobs(r)
    assert j["stranded_signal"] == 1 and "bigwig" not in j
    assert "stranded_signal" not in jobs(dry_run(ANA, outputDir=str(out), stranded_plots="False"))


@needs_snakemake
def test_stranded_config_values_reach_the_command_and_bad_formats_are_rejected():
    r = dry_run(FULL, track_loci="[Gapdh, 'chr1:100-200']", track_flank="2000", track_bin_size="50", track_max_loci="7",
                metagene_flank="1500", track_static_formats="[svg]")
    assert r.returncode == 0, r.stderr
    r2 = subprocess.run([SNAKEMAKE, "-n", "-p", "-s", "workflows/TTpipe.snakefile", "--cores", "1", "--configfile", MAIN, FULL,
                         "--config", "track_loci=[Gapdh, 'chr1:100-200']", "track_flank=2000", "track_bin_size=50",
                         "track_max_loci=7", "metagene_flank=1500", "track_static_formats=[svg]"],
                        cwd=ROOT, capture_output=True, text=True)
    cmd = r2.stdout.split("stranded_signal.py")[1].split("2>")[0]
    for frag in ("--loci Gapdh chr1:100-200", "--flank 2000", "--bin-size 50", "--max-loci 7", "--metagene-flank 1500",
                 "--static-formats svg", "--top-degs 10", "--degs "):
        assert frag in cmd, frag
    bad = dry_run(FULL, track_static_formats="[gif]")
    assert bad.returncode != 0 and "track_static_formats" in (bad.stderr + bad.stdout)
