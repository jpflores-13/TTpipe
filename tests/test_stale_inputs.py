"""Editing the samplesheet must regenerate everything that is derived from it.

Regression test for a bug seen on a real run: the sample table (samples_deseq2.txt) did not list the samplesheet
as an input, so after samples were added Snakemake kept the old table and DESeq2 failed with
'sample names differ between the count matrix and the sample table'.

The `sample_table` rule is pure Python (no external tools), so this test really runs it.
"""
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SNAKEMAKE = shutil.which("snakemake") or str(Path(sys.executable).parent / "snakemake")
pytestmark = pytest.mark.skipif(not Path(SNAKEMAKE).exists(), reason="snakemake not installed")


def run(args, outdir, sheet):
    cmd = [SNAKEMAKE, "-s", "workflows/TTpipe.snakefile", "--cores", "1",
           *args,                  # targets / flags go BEFORE --configfile and --config (both swallow following values)
           "--configfile", "config/TTpipeConfig.yaml", "tests/dryrun/config_analysis_only.yaml",
           "--config", f"samplesheet={sheet}", f"outputDir={outdir}"]
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)


def test_sample_table_follows_the_samplesheet(tmp_path):
    sheet = tmp_path / "samplesheet.txt"
    shutil.copy(ROOT / "tests/dryrun/samplesheet.txt", sheet)
    out = tmp_path / "out"
    target = out / "counts" / "samples_deseq2.txt"

    r = run([str(target)], out, sheet)
    assert r.returncode == 0, r.stderr
    assert len(target.read_text().splitlines()) == 1 + 4          # header + 4 samples

    # unchanged samplesheet: nothing to redo
    r = run(["-n", str(target)], out, sheet)
    assert "Nothing to be done" in (r.stdout + r.stderr)

    # add two samples (a new condition); wait so the samplesheet's mtime is newer than the table's
    time.sleep(1.5)
    with open(sheet, "a") as fh:
        fh.write("CI_6H_rep1\t6H\trep1\tR1.fastq.gz\tR2.fastq.gz\ttests/dryrun/\n")
        fh.write("CI_6H_rep2\t6H\trep2\tR1.fastq.gz\tR2.fastq.gz\ttests/dryrun/\n")

    r = run([str(target)], out, sheet)
    assert r.returncode == 0, r.stderr
    rows = target.read_text().splitlines()
    assert len(rows) == 1 + 6
    assert any(line.startswith("CI_6H_rep2\t6H\trep2") for line in rows)


def test_report_rules_list_the_samplesheet_as_input():
    cmd = [SNAKEMAKE, "-n", "-s", "workflows/TTpipe.snakefile", "--cores", "1",
           "--configfile", "config/TTpipeConfig.yaml", "tests/dryrun/config.yaml"]
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    blocks = r.stdout.split("\n\n")
    for rule in ("sample_table", "qc_summary", "qc_report"):
        block = next((b for b in blocks if f"rule {rule}:" in b), "")
        assert "tests/dryrun/samplesheet.txt" in block.split("output:")[0], rule
