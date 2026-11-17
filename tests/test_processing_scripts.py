"""Unit tests for the processing helper scripts (tiny synthetic inputs)."""
import gzip
import io
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import gene_info as gi  # noqa: E402
import make_count_matrix as mcm  # noqa: E402
import scale_factors as sf  # noqa: E402

GENCODE = (
    '#comment\n'
    'chr1\tHAVANA\tgene\t100\t199\t.\t+\t.\tgene_id "G1.5"; gene_type "protein_coding"; gene_name "Aa";\n'
    'chr1\tHAVANA\ttranscript\t100\t199\t.\t+\t.\tgene_id "G1.5"; transcript_id "T1";\n'
    'chr2\tHAVANA\tgene\t10\t10\t.\t-\t.\tgene_id "G2.1"; gene_type "lncRNA"; gene_name "Bb";\n'
)
ENSEMBL = (
    '1\tensembl\tgene\t5\t14\t.\t-\t.\tgene_id "ENSG1"; gene_version "1"; gene_name "Cc"; gene_biotype "miRNA";\n'
    '1\tensembl\tgene\t50\t60\t.\t+\t.\tgene_id "ENSG2"; gene_biotype "rRNA";\n'
)


def test_gene_info_gencode():
    rows = gi.read_gene_info(GENCODE.splitlines(True))
    assert [r["gene_id"] for r in rows] == ["G1.5", "G2.1"]
    assert rows[0] == {"gene_id": "G1.5", "gene_name": "Aa", "gene_type": "protein_coding", "chr": "chr1",
                       "start": 100, "end": 199, "strand": "+", "length": 100}
    assert rows[1]["length"] == 1 and rows[1]["strand"] == "-"


def test_gene_info_ensembl_and_missing_name():
    rows = gi.read_gene_info(ENSEMBL.splitlines(True))
    assert rows[0]["gene_type"] == "miRNA" and rows[0]["gene_name"] == "Cc"
    assert rows[1]["gene_name"] == "ENSG2" and rows[1]["gene_type"] == "rRNA"


def test_gene_info_cli_gzip(tmp_path):
    gtf = tmp_path / "a.gtf.gz"
    with gzip.open(gtf, "wt") as fh:
        fh.write(GENCODE)
    out = tmp_path / "gi.tsv"
    gi.main(["--gtf", str(gtf), "--out", str(out)])
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == gi.COLUMNS and len(df) == 2


def _htseq(path, d):
    path.write_text("".join(f"{k}\t{v}\n" for k, v in d.items()))
    return str(path)


@pytest.fixture
def htseq_files(tmp_path):
    summ = {"__no_feature": 5, "__ambiguous": 1, "__too_low_aQual": 0, "__not_aligned": 0, "__alignment_not_unique": 2}
    f1 = _htseq(tmp_path / "f1", {"gP": 10, "gM": 1, "gN": 3, **summ})
    f2 = _htseq(tmp_path / "f2", {"gP": 20, "gM": 2, "gN": 4, **summ})
    r1 = _htseq(tmp_path / "r1", {"gP": 1, "gM": 7, "gN": 5, **summ})
    r2 = _htseq(tmp_path / "r2", {"gP": 2, "gM": 8, "gN": 6, **summ})
    gene_info = tmp_path / "gene_info.tsv"
    gene_info.write_text("gene_id\tstrand\ngP\t+\ngM\t-\ngN\t.\n")
    return [f1, f2], [r1, r2], str(gene_info), tmp_path


def test_matrix_drops_summary_rows(htseq_files):
    f, r, _, _ = htseq_files
    m = mcm.build_matrix(f, ["s1", "s2"])
    assert list(m.index) == ["gP", "gM", "gN"] and list(m.columns) == ["s1", "s2"]
    assert m.loc["gP", "s2"] == 20


def test_sum_and_sense(htseq_files):
    f, r, g, _ = htseq_files
    fwd, rev = mcm.build_matrix(f, ["s1", "s2"]), mcm.build_matrix(r, ["s1", "s2"])
    s = mcm.deseq2_matrix(fwd, rev, "sum")
    assert s.loc["gP"].tolist() == [11, 22]
    sense = mcm.deseq2_matrix(fwd, rev, "sense", mcm.read_gene_info(g))
    assert sense.loc["gP"].tolist() == [10, 20]      # + gene: forward
    assert sense.loc["gM"].tolist() == [7, 8]        # - gene: reverse
    assert sense.loc["gN"].tolist() == [8, 10]       # unknown strand: sum
    assert all(dt.kind == "i" for dt in sense.dtypes)


def test_sense_missing_gene_info(htseq_files):
    f, r, _, tmp = htseq_files
    fwd, rev = mcm.build_matrix(f, ["s1", "s2"]), mcm.build_matrix(r, ["s1", "s2"])
    gi_df = pd.DataFrame({"strand": ["+"]}, index=pd.Index(["gP"], name="gene_id"))
    with pytest.raises(ValueError):
        mcm.deseq2_matrix(fwd, rev, "sense", gi_df)


def test_cli_outputs(htseq_files):
    f, r, g, tmp = htseq_files
    mcm.main(["--names", "s1", "s2", "--forward", *f, "--reverse", *r, "--count-strand", "sense",
              "--gene-info", g, "--out-forward", str(tmp / "F"), "--out-reverse", str(tmp / "R"),
              "--out-deseq2", str(tmp / "D")])
    lines = (tmp / "D").read_text().splitlines()
    assert lines[0] == "gene_id\ts1\ts2" and lines[1] == "gP\t10\t20"
    assert (tmp / "F").read_text().splitlines()[0] == "gene_id\ts1\ts2"


def test_gene_set_mismatch(tmp_path):
    a = _htseq(tmp_path / "a", {"g1": 1})
    b = _htseq(tmp_path / "b", {"g2": 1})
    with pytest.raises(ValueError):
        mcm.build_matrix([a, b], ["a", "b"])


def test_scale_factors_computed():
    rows = sf.compute_scale_factors({"DMSO": 800, "2H": 1000, "6H": 500})
    d = {r[0]: r for r in rows}
    assert d["2H"][2] == 1 and d["DMSO"][2] == pytest.approx(1.25) and d["6H"][2] == 2
    assert {r[3] for r in rows} == {"computed"}


def test_scale_factors_override_and_errors():
    rows = sf.compute_scale_factors({"A": 100, "B": 200}, {"A": 1.5})
    assert rows[0][2] == 1.5 and rows[0][3] == "config" and rows[1][2] == 1 and rows[1][3] == "computed"
    with pytest.raises(ValueError):
        sf.compute_scale_factors({"A": 100}, {"Z": 2})
    with pytest.raises(ValueError):
        sf.compute_scale_factors({"A": 0, "B": 5})


def test_scale_factors_cli(tmp_path):
    p1, p2 = tmp_path / "a.txt", tmp_path / "b.txt"
    p1.write_text("S1\t100\n")
    p2.write_text("S2\t400\n")
    out = tmp_path / "sf.tsv"
    sf.main(["--pairs", str(p1), str(p2), "--out", str(out), "--override", "S1=2"])
    lines = out.read_text().splitlines()
    assert lines[0] == "condition\tmapped_pairs\tscale_factor\tsource"
    assert lines[1] == "S1\t100\t2\tconfig" and lines[2] == "S2\t400\t1\tcomputed"
