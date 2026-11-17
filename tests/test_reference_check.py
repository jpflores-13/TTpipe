"""The chromosome-naming check in common.smk (the guide's Ensembl-FASTA + GENCODE-GTF trap)."""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SNAKEMAKE = shutil.which("snakemake") or str(Path(sys.executable).parent / "snakemake")
pytestmark = pytest.mark.skipif(not Path(SNAKEMAKE).exists(), reason="snakemake not installed")


def dry_run(*extra):
    cmd = [SNAKEMAKE, "-n", "-s", "workflows/TTpipe.snakefile", "--cores", "1",
           "--configfile", "config/TTpipeConfig.yaml", "tests/dryrun/config.yaml", *extra]
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)


def test_mismatched_index_names_are_rejected(tmp_path):
    idx = tmp_path / "idx"
    idx.mkdir()
    (idx / "chrNameLength.txt").write_text("1\t1000\n2\t900\n")      # Ensembl style; fixture GTF uses 'chrP'
    r = dry_run("--config", f"star_index={idx}")
    assert r.returncode != 0
    assert "chromosome names differ" in (r.stderr + r.stdout)


def test_mismatched_fasta_is_rejected_when_building_index(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">1 dna:chromosome\nACGT\n")
    r = dry_run("--config", f"genome_fasta={fa}", f"star_index={tmp_path / 'nonexistent'}", "build_star_index=True")
    assert r.returncode != 0
    assert "chromosome names differ" in (r.stderr + r.stdout)


def test_matching_names_pass():
    r = dry_run()
    assert r.returncode == 0, r.stderr
