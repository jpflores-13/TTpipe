"""Optional spike-in module (EXPERIMENTAL): scripts on tiny hand-made inputs, and Snakemake dry runs.

The spike-in fixtures in tests/dryrun/ (spike.fa, spike.gtf) mirror the real ERCC92 files: one contig per spike-in,
a GTF with `exon` lines only, gene_id equal to the contig name and a separate transcript accession. The real ERCC
files are not in the repository.
"""
import csv
import subprocess
import sys
from pathlib import Path

import pytest

from test_dryrun import FULL, dry_run, jobs, needs_snakemake

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
SPIKE_FA = "tests/dryrun/spike.fa"
SPIKE_GTF = "tests/dryrun/spike.gtf"


def run(script, *args):
    return subprocess.run([sys.executable, str(SCRIPTS / script), *map(str, args)], capture_output=True, text=True)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------- build_spike_reference.py ----------------

def test_combined_reference_prefixes_and_synthesises_gene_lines(tmp_path):
    out = {k: tmp_path / k for k in ("fa", "gtf", "contigs", "genes")}
    r = run("build_spike_reference.py", "--host-fasta", ROOT / "tests/dryrun/genome.fa",
            "--host-gtf", ROOT / "tests/dryrun/genes.gtf", "--spike-fasta", ROOT / SPIKE_FA,
            "--spike-gtf", ROOT / SPIKE_GTF, "--prefix", "spike_", "--out-fasta", out["fa"], "--out-gtf", out["gtf"],
            "--out-contigs", out["contigs"], "--out-genes", out["genes"])
    assert r.returncode == 0, r.stderr
    headers = [ln[1:].strip() for ln in out["fa"].read_text().splitlines() if ln.startswith(">")]
    assert headers == ["chrP", "spike_ERCC-00002", "spike_ERCC-00003"]
    gtf = [ln.split("\t") for ln in out["gtf"].read_text().splitlines() if ln and not ln.startswith("#")]
    spike = [f for f in gtf if f[0].startswith("spike_")]
    genes = [f for f in spike if f[2] == "gene"]
    assert len(genes) == 2 and all('gene_id "spike_ERCC-' in f[8] for f in genes)    # HTSeq counts --type gene
    assert all('transcript_id "spike_DQ' in f[8] for f in spike if f[2] == "exon")   # transcript ids prefixed too
    assert [r["contig"] for r in read_tsv(out["contigs"])] == ["spike_ERCC-00002", "spike_ERCC-00003"]
    assert {r["gene_id"] for r in read_tsv(out["genes"])} == {"spike_ERCC-00002", "spike_ERCC-00003"}
    host = [f for f in gtf if f[0] == "chrP"]
    assert len(host) == 2                                                           # host GTF unchanged


def test_host_contig_with_the_prefix_is_refused(tmp_path):
    fa = tmp_path / "host.fa"
    fa.write_text(">spike_chr1\nACGT\n")
    r = run("build_spike_reference.py", "--host-fasta", fa, "--host-gtf", ROOT / "tests/dryrun/genes.gtf",
            "--spike-fasta", ROOT / SPIKE_FA, "--prefix", "spike_", "--out-fasta", tmp_path / "o.fa",
            "--out-gtf", tmp_path / "o.gtf", "--out-contigs", tmp_path / "o.tsv")
    assert r.returncode != 0 and "spike_" in r.stderr


# ---------------- spike_summary.py -> scale_factors.py / spike_size_factors.py ----------------

def idxstats(path, host, ercc2, ercc3):
    path.write_text(f"chrP\t80\t{host}\t0\nspike_ERCC-00002\t40\t{ercc2}\t0\nspike_ERCC-00003\t40\t{ercc3}\t0\n*\t0\t0\t5\n")
    return path


@pytest.fixture
def spike_counts(tmp_path):
    # records (not pairs): A1 has twice the spike-in of the others
    a1 = idxstats(tmp_path / "a1.txt", 2000, 400, 0)
    a2 = idxstats(tmp_path / "a2.txt", 2000, 200, 0)
    b1 = idxstats(tmp_path / "b1.txt", 4000, 200, 0)
    b2 = idxstats(tmp_path / "b2.txt", 4000, 200, 0)
    out, per = tmp_path / "spike_counts.tsv", tmp_path / "contigs.tsv"
    r = run("spike_summary.py", "--sample", "A1", "DMSO", a1, "--sample", "A2", "DMSO", a2,
            "--sample", "B1", "2H", b1, "--sample", "B2", "2H", b2, "--prefix", "spike_",
            "--out-counts", out, "--out-contigs", per)
    assert r.returncode == 0, r.stderr
    return out


def test_spike_summary_counts_pairs_per_group(spike_counts):
    rows = {r["sample"]: r for r in read_tsv(spike_counts)}
    assert rows["A1"]["host_pairs"] == "1000" and rows["A1"]["spike_pairs_all"] == "200"
    assert rows["B1"]["host_pairs"] == "2000" and rows["B1"]["spike_pairs_used"] == "100"
    assert float(rows["A1"]["spike_pct"]) == pytest.approx(200 / 1200 * 100, abs=1e-3)   # written with 4 decimals


def test_condition_track_factors_from_spike_and_from_host(spike_counts, tmp_path):
    out = tmp_path / "sf.tsv"
    r = run("scale_factors.py", "--spike-counts", spike_counts, "--group-by", "condition", "--basis", "spike",
            "--order", "DMSO", "2H", "--out", out)
    assert r.returncode == 0, r.stderr
    sf = {r["condition"]: r for r in read_tsv(out)}
    assert float(sf["DMSO"]["scale_factor"]) == 1 and float(sf["2H"]["scale_factor"]) == pytest.approx(1.5)
    assert sf["2H"]["source"] == "spike-in"                          # DMSO 300 spike pairs, 2H 200 -> 300 / 200
    r = run("scale_factors.py", "--spike-counts", spike_counts, "--group-by", "condition", "--basis", "host",
            "--order", "DMSO", "2H", "--out", out)
    sf = {r["condition"]: r for r in read_tsv(out)}
    assert float(sf["DMSO"]["scale_factor"]) == pytest.approx(2) and float(sf["2H"]["scale_factor"]) == 1


def test_sample_track_factor(spike_counts):
    r = run("scale_factors.py", "--spike-counts", spike_counts, "--group-by", "sample", "--basis", "spike",
            "--factor-of", "B1")
    assert r.returncode == 0, r.stderr
    assert float(r.stdout) == pytest.approx(2)                       # A1 200 pairs / B1 100


def test_deseq2_total_size_factors(spike_counts, tmp_path):
    out = tmp_path / "size_factors.tsv"
    r = run("spike_size_factors.py", "--method", "total", "--spike-counts", spike_counts, "--out", out)
    assert r.returncode == 0, r.stderr
    sf = {r["sample"]: float(r["size_factor"]) for r in read_tsv(out)}
    gm = (200 * 100 * 100 * 100) ** 0.25
    assert sf["A1"] == pytest.approx(200 / gm) and sf["B2"] == pytest.approx(100 / gm)


def test_zero_spike_in_pairs_is_a_clear_error(tmp_path):
    a = idxstats(tmp_path / "a.txt", 2000, 0, 0)
    b = idxstats(tmp_path / "b.txt", 2000, 200, 0)
    counts = tmp_path / "c.tsv"
    run("spike_summary.py", "--sample", "A", "DMSO", a, "--sample", "B", "2H", b, "--prefix", "spike_",
        "--out-counts", counts, "--out-contigs", tmp_path / "p.tsv")
    r = run("spike_size_factors.py", "--method", "total", "--spike-counts", counts, "--out", tmp_path / "o.tsv")
    assert r.returncode != 0 and "zero spike-in" in r.stderr


# ---------------- dry runs ----------------

def spike_config(tmp_path, **spike):
    block = {"enabled": True, "fasta": SPIKE_FA, "gtf": SPIKE_GTF, **spike}
    lines = ["build_star_index: true", "star_index: 'output_ci/reference/star_index_spike'", "spike_in:"]
    for k, v in block.items():
        val = ("true" if v else "false") if isinstance(v, bool) else (repr(v) if isinstance(v, list) else f"'{v}'")
        lines.append(f"  {k}: {val}")
    p = tmp_path / "spike.yaml"
    p.write_text("\n".join(lines) + "\n")
    return str(p)


@needs_snakemake
def test_disabled_is_the_plain_workflow():
    j = jobs(dry_run(FULL))
    assert not any(k.startswith("spike") for k in j)


@needs_snakemake
def test_enabled_total_method(tmp_path):
    r = dry_run(FULL, spike_config(tmp_path))
    assert r.returncode == 0, r.stderr
    j = jobs(r)
    assert j["spike_reference"] == 1 and j["star_index"] == 1
    assert j["spike_idxstats"] == 4 and j["spike_summary"] == 1 and j["spike_size_factors"] == 1
    assert j["scale_factors"] == 1 and j["count_matrix"] == 1 and j["deseq2"] == 1
    assert "spike/size_factors.tsv" in r.stdout and "counts/counts_spike.txt" in r.stdout
    assert "deseq2_size_factors.tsv" in r.stdout


@needs_snakemake
def test_enabled_genes_method_uses_control_genes(tmp_path):
    r = dry_run(FULL, spike_config(tmp_path, size_factor_method="genes"))
    assert r.returncode == 0, r.stderr
    assert "spike/control_genes.txt" in r.stdout and "spike/size_factors.tsv" not in r.stdout


@needs_snakemake
def test_count_only_has_no_size_factor_job(tmp_path):
    r = dry_run(FULL, spike_config(tmp_path, normalize=[]))
    assert r.returncode == 0, r.stderr
    j = jobs(r)
    assert j["spike_summary"] == 1 and "spike_size_factors" not in j
    assert "deseq2_size_factors.tsv" not in r.stdout


@needs_snakemake
def test_existing_index_without_spike_contigs_is_refused(tmp_path):
    cfg = spike_config(tmp_path)
    Path(cfg).write_text(Path(cfg).read_text().replace("build_star_index: true", "build_star_index: false")
                         .replace("output_ci/reference/star_index_spike", "tests/dryrun/star_index"))
    r = dry_run(FULL, cfg)
    assert r.returncode != 0 and "does not contain the spike-in contigs" in (r.stderr + r.stdout)
