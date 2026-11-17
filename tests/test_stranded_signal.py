"""Tests for scripts/stranded_signal.py on SYNTHETIC data built in this file.

The genome, genes, bigWigs and counts below are INVENTED (one chromosome 'chrS', a + gene and a - gene with
hand-chosen signal values, plus filler genes so the metagene groups have enough members). Passing these tests shows
that strand orientation, locus resolution, binning, gene-model parsing, the strand-table arithmetic and the page
JavaScript behave as designed on those inputs; it says nothing about real TT-seq data. Tests that need pyBigWig are
skipped when it is not importable; the browser test runs in headless Chrome when Chrome is installed.

Synthetic signal recipe (per condition, multiplied by k = 1 for condition A and 2 for condition B)
  a gene on strand s with TSS at its 5' end:
    sense strand bigWig : 10 over the gene body, 40 over the first 100 bp after the TSS (a spike)
    antisense bigWig    : 1 over the gene body and 7 over the 500 bp UPSTREAM of the TSS (a bump)
  so sense and antisense have different, known shapes, and a - gene's numbers sit in the reverse bigWig.
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

import stranded_signal as ss  # noqa: E402
from gene_info import read_gene_info, write_gene_info  # noqa: E402

pyBigWig = pytest.importorskip("pyBigWig", reason="pyBigWig not installed: stranded-signal data tests skipped")

CHROM, CLEN = "chrS", 600_000
K = {"A": 1, "B": 2}
CONDS = ["A", "B"]
SAMPLES = [("A_1", "A"), ("A_2", "A"), ("B_1", "B"), ("B_2", "B")]
# (gene_id, name, type, strand, start (1-based), end)
GENES = [("SYNP0001.3", "SynP", "protein_coding", "+", 30001, 38000),
         ("SYNM0001.5", "SynM", "lncRNA", "-", 60001, 70000)]
for i in range(25):
    s = 100_001 + i * 10_000
    GENES.append((f"PC{i:02d}", f"Pc{i}", "protein_coding", "+" if i % 2 == 0 else "-", s, s + 1999))
for i in range(22):
    s = 350_001 + i * 10_000
    GENES.append((f"LN{i:02d}", f"Ln{i}", "lncRNA", "+" if i % 2 == 0 else "-", s, s + 1999))
LOW_EXPR = {"LN20", "LN21"}        # two lncRNA fillers get (almost) no counts: not expressed


def gtf_text():
    lines = []
    for gid, name, typ, strand, s, e in GENES:
        attr = f'gene_id "{gid}"; gene_name "{name}"; gene_type "{typ}";'
        lines.append(f"{CHROM}\tsyn\tgene\t{s}\t{e}\t.\t{strand}\t.\t{attr}")
        if gid == "SYNP0001.3":
            exons = [(30001, 30500), (33001, 33800), (37001, 38000)]
        elif gid == "SYNM0001.5":
            exons = [(60001, 61000), (64001, 64800), (69001, 70000)]
        else:
            exons = [(s, e)]
        for xs, xe in exons:
            lines.append(f"{CHROM}\tsyn\texon\t{xs}\t{xe}\t.\t{strand}\t.\t{attr} transcript_id \"T{gid}\";")
        lines.append(f"{CHROM}\tsyn\texon\t{s}\t{s + 5}\t.\t{strand}\t.\t{attr} transcript_id \"T2{gid}\";")  # second transcript, overlapping exon
    lines.append(f"chrOther\tsyn\texon\t1\t100\t.\t+\t.\tgene_id \"FAR\"; gene_name \"Far\";")
    return "\n".join(lines) + "\n"


def intervals(strand_file):
    """Sorted (start0, end0, value) triples (k = 1) of the forward or reverse bigWig for ALL genes."""
    out = []
    for gid, name, typ, strand, s, e in GENES:
        s0, e0 = s - 1, e
        sense_file = "forward" if strand == "+" else "reverse"
        if strand_file == sense_file:     # the gene's own strand
            if strand == "+":
                out += [(s0, s0 + 100, 40.0), (s0 + 100, e0, 10.0)]
            else:
                out += [(s0, e0 - 100, 10.0), (e0 - 100, e0, 40.0)]
        else:                             # the opposite strand: body 1 and a 500 bp bump upstream of the TSS
            if strand == "+":
                out += [(s0 - 500, s0, 7.0), (s0, e0, 1.0)]
            else:
                out += [(s0, e0, 1.0), (e0, e0 + 500, 7.0)]
    return sorted(out)


def write_bw(path, items, k, nan_gap=False):
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([(CHROM, CLEN), ("chrOther", 1000)])
    items = sorted(items)
    bw.addEntries([CHROM] * len(items), [a for a, _, _ in items], ends=[b for _, b, _ in items],
                  values=[v * k for _, _, v in items])
    bw.close()


@pytest.fixture(scope="module")
def syn(tmp_path_factory):
    root = tmp_path_factory.mktemp("stranded")
    (root / "tracks").mkdir()
    (root / "counts").mkdir()
    bws = {}
    for c in CONDS:
        bws[c] = {}
        for strand in ("forward", "reverse"):
            p = root / "tracks" / f"{c}_{strand}.bw"
            write_bw(p, intervals(strand), K[c])
            bws[c][strand] = str(p)
    gtf = root / "genes.gtf"
    gtf.write_text(gtf_text())
    with open(gtf) as fh:
        rows = read_gene_info(fh)
    gi_path = root / "counts" / "gene_info.tsv"
    with open(gi_path, "w") as fh:
        write_gene_info(rows, fh)
    # counts: sense 100*(i+1) per sample, antisense 10 (so the antisense fraction is 10/(100(i+1)+10) per gene)
    ids = [g[0] for g in GENES]
    f = pd.DataFrame(0, index=ids, columns=[s for s, _ in SAMPLES], dtype=float)
    r = f.copy()
    for i, (gid, name, typ, strand, s, e) in enumerate(GENES):
        sense = 100.0 * (i + 1) if gid not in LOW_EXPR else 1.0
        anti = 10.0 if gid not in LOW_EXPR else 0.0
        for smp, _ in SAMPLES:
            (f if strand == "+" else r).loc[gid, smp] = sense
            (r if strand == "+" else f).loc[gid, smp] = anti
    for name, m in (("counts_forward.txt", f), ("counts_reverse.txt", r), ("counts_deseq2.txt", f + r)):
        m.index.name = "gene_id"
        m.to_csv(root / "counts" / name, sep="\t")
    pd.DataFrame({"sample": [s for s, _ in SAMPLES], "condition": [c for _, c in SAMPLES], "replicate": "r"}).to_csv(
        root / "counts" / "samples_deseq2.txt", sep="\t", index=False)
    pd.DataFrame({"condition": CONDS, "mapped_pairs": [1000, 500], "scale_factor": [1.0, 2.0], "source": "computed"}).to_csv(
        root / "tracks" / "scale_factors.tsv", sep="\t", index=False)
    # DESeq2-like tables: SYNM most significant in contrast 1, a NA-padj gene in contrast 2
    d = root / "deseq2"
    d.mkdir()
    for name, order in (("B_vs_A", ["PC03", "SYNM0001.5", "SYNP0001.3", "PC01"]), ("A_vs_B", ["SYNP0001.3", "PC01"])):
        padj = {g: 1e-8 * 10.0 ** k for k, g in enumerate(order)}   # first = most significant
        df = pd.DataFrame({"gene_id": order + ["PC05", "PC07"], "baseMean": 5.0, "log2FoldChange": 1.0, "lfcSE": .1,
                           "stat": 1.0, "pvalue": 0.5, "padj": [padj[g] for g in order] + [np.nan, 0.9]})
        df.to_csv(d / f"degs_{name}.csv", index=False)
    return {"root": root, "bws": bws, "gtf": str(gtf), "gi": str(gi_path), "gene_table": ss.read_gene_table(str(gi_path)),
            "counts_f": str(root / "counts" / "counts_forward.txt"), "counts_r": str(root / "counts" / "counts_reverse.txt"),
            "counts_m": str(root / "counts" / "counts_deseq2.txt"), "samples": str(root / "counts" / "samples_deseq2.txt"),
            "sf": str(root / "tracks" / "scale_factors.tsv"),
            "degs": {"B_vs_A": str(d / "degs_B_vs_A.csv"), "A_vs_B": str(d / "degs_A_vs_B.csv")}}


def prep(syn, **kw):
    args = dict(gtf=syn["gtf"], scale_factors=syn["sf"], counts_forward=syn["counts_f"], counts_reverse=syn["counts_r"],
                counts_matrix=syn["counts_m"], samples=syn["samples"], banner="SYNTHETIC TEST DATA")
    args.update(kw)
    return ss.prepare_data(syn["bws"], CONDS, syn["gene_table"], **args)


# ---------------- gene table and coordinate convention ----------------

def test_gene_table_keeps_gtf_coordinates(syn):
    gt = syn["gene_table"]
    row = gt[gt["gene_id"] == "SYNP0001.3"].iloc[0]
    assert (row["start"], row["end"], row["strand"]) == (30001, 38000, "+")      # 1-based inclusive, as gene_info.py writes them
    loc = ss.gene_locus(row, "t")
    assert (loc["start0"], loc["end0"]) == (30000, 38000)                          # bigWig coordinates: 0-based half-open


# ---------------- locus resolution ----------------

def test_resolve_symbol_case_insensitive_id_with_and_without_version(syn):
    gt = syn["gene_table"]
    loci, nf, notes = ss.select_loci(["synp", "SYNP0001", "SYNM0001.5", "SynM"], gt)
    assert nf == []
    assert [l["gene_id"] for l in loci] == ["SYNP0001.3", "SYNM0001.5"]          # de-duplicated, order kept
    assert loci[0]["sources"] == ["config: synp", "config: SYNP0001"]            # both ways of asking are remembered
    assert loci[0]["label"] == "SynP"


def test_resolve_id_without_version_and_exact_with_other_version(syn):
    gt = syn["gene_table"]
    loci, nf, _ = ss.select_loci(["SYNP0001"], gt)
    assert [l["gene_id"] for l in loci] == ["SYNP0001.3"]
    loci, nf, _ = ss.select_loci(["SYNP0001.99"], gt)       # a different version number still resolves by the unversioned id
    assert [l["gene_id"] for l in loci] == ["SYNP0001.3"]


def test_ambiguous_symbol_takes_all_and_says_so(syn):
    gt = pd.concat([syn["gene_table"], syn["gene_table"].iloc[[0]].assign(gene_id="SYNP9999.1", start=500, end=900)],
                   ignore_index=True)
    loci, nf, notes = ss.select_loci(["SynP"], gt)
    assert {l["gene_id"] for l in loci} == {"SYNP0001.3", "SYNP9999.1"}
    assert len(notes) == 1 and "matches 2 genes" in notes[0] and "SYNP9999.1" in notes[0]
    assert all("(" in l["label"] for l in loci)                                    # labels disambiguated by gene_id


def test_resolve_coordinates_and_not_found_is_not_fatal(syn):
    gt = syn["gene_table"]
    loci, nf, _ = ss.select_loci(["chrS:1,000-2,000", "NoSuchGene", "chrS:50-10", "chrS:x"], gt)
    assert len(loci) == 1 and (loci[0]["chrom"], loci[0]["start0"], loci[0]["end0"]) == ("chrS", 999, 2000)
    assert loci[0]["label"] == "chrS:1000-2000"
    assert [n["entry"] for n in nf] == ["NoSuchGene", "chrS:50-10", "chrS:x"]       # reported, in order, nothing raised
    assert all(n["reason"] for n in nf)


def test_resolve_bed_with_and_without_names(tmp_path, syn):
    bed = tmp_path / "r.bed"
    bed.write_text("track name=x\n# comment\nchrS\t1000\t1500\tenhA\t0\t+\nchrS\t3000\t3400\nchrS\t5\t5\nbad line\n")
    loci, nf, _ = ss.select_loci([], syn["gene_table"], bed=str(bed))
    assert [l["label"] for l in loci] == ["enhA", "chrS:3001-3400"]
    assert loci[0]["strand"] == "+" and loci[1]["strand"] == ""
    assert (loci[0]["start0"], loci[0]["end0"]) == (1000, 1500)                     # BED is already 0-based half-open
    assert ss.select_loci([], syn["gene_table"], bed=str(tmp_path / "missing.bed"))[1][0]["reason"].startswith("track_regions_bed")


def test_top_degs_per_contrast_smallest_padj_and_na_skipped(syn):
    loci, nf, _ = ss.select_loci([], syn["gene_table"], deg_tables=syn["degs"], top_n=2)
    ids = [l["gene_id"] for l in loci]
    assert ids == ["PC03", "SYNM0001.5", "SYNP0001.3", "PC01"]                      # B_vs_A: PC03, SYNM; A_vs_B: SYNP, PC01
    assert ss.top_deg_ids(syn["degs"]["B_vs_A"], 10) == ["PC03", "SYNM0001.5", "SYNP0001.3", "PC01", "PC07"]   # PC05 has NA padj
    assert nf == []
    assert ss.select_loci([], syn["gene_table"], deg_tables=syn["degs"], top_n=0)[0] == []
    assert ss.select_loci([], syn["gene_table"], deg_tables={"x": "/no/such.csv"}, top_n=3)[0] == []   # missing table: silently skipped


def test_combined_dedup_and_cap(tmp_path, syn):
    bed = tmp_path / "r.bed"
    bed.write_text("chrS\t1000\t1500\tenhA\nchrS\t3000\t3400\tenhB\n")
    loci, nf, notes = ss.select_loci(["SynP", "chrS:1001-1500", "Nope"], syn["gene_table"], deg_tables=syn["degs"],
                                     top_n=2, bed=str(bed), max_loci=4)
    assert len(loci) == 4 and [l["label"] for l in loci][:2] == ["SynP", "chrS:1001-1500"]
    # chrS:1001-1500 equals the BED region enhA (1000-1500 0-based) and is not drawn twice
    full = ss.select_loci(["SynP", "chrS:1001-1500", "Nope"], syn["gene_table"], deg_tables=syn["degs"], top_n=2, bed=str(bed))[0]
    assert len(full) == 2 + 3 + 1                                                    # SynP (also a top DEG: merged), coord/enhA merged, PC03, SynM, PC01, enhB
    assert [n["entry"] for n in nf] == ["Nope"]
    assert any("only the first 4" in n for n in notes)
    assert "top DEG: A_vs_B" in full[0]["sources"]                                  # SynP was asked for AND is a top DEG: one locus, both sources


# ---------------- bins and NaN ----------------

def test_plan_bins_whole_bins_and_cap():
    assert ss.plan_bins(1000, 2050, 10**6, 100) == (1000, 11, 100)
    s, n, bs = ss.plan_bins(0, 10_000_000, 10**8, 100)
    assert bs % 100 == 0 and n <= ss.MAX_BINS and n * bs >= 10_000_000 - bs
    assert ss.plan_bins(999_950, 1_000_300, 1_000_000, 100)[1:] == (1, 50)
    s, n, bs = ss.plan_bins(999_950, 1_000_300, 1_000_000, 100)
    assert s + n * bs <= 1_000_000                                                   # never past the chromosome end


def test_bin_mean_counts_uncovered_bases_as_zero(tmp_path):
    p = tmp_path / "gap.bw"
    bw = pyBigWig.open(str(p), "w")
    bw.addHeader([(CHROM, 10_000)])
    bw.addEntries([CHROM] * 3, [0, 100, 300], ends=[50, 200, 400], values=[2.0, 4.0, 6.0])
    bw.close()
    bw = pyBigWig.open(str(p))
    got = ss.bin_mean(bw, CHROM, 0, 500, 5)
    assert got.shape == (5,)
    # bin 0: 50 bp of value 2 in 100 bp -> 1.0 (pyBigWig's own `mean` would say 2.0); bin 2: nothing -> 0 (None -> 0)
    assert got[0] == pytest.approx(1.0) and got[1] == pytest.approx(4.0) and got[2] == 0.0
    assert got[3] == pytest.approx(6.0) and got[4] == 0.0
    assert bw.stats(CHROM, 0, 100, type="mean", nBins=1)[0] == 2.0                  # documents why `sum` / width is used
    assert not np.isnan(got).any()
    bw.close()


# ---------------- orientation: locus tracks ----------------

def test_plus_gene_sense_is_forward_minus_gene_sense_is_reverse(syn):
    data = prep(syn, loci_entries=["SynP", "SynM"], flank=5000, bin_size=100)
    p, m = data["loci"]
    assert (p["label"], m["label"]) == ("SynP", "SynM")
    # SynP (+, 30001-38000): window 25000-43000 (0-based 24999 start), bins of 100 -> index of gene start:
    # window start0 = 25000 - 1 = 24999?  start0 = 30000 - 5000 = 25000, so the gene starts at bin 50
    assert p["win_start"] == 25001 and p["bin"] == 100 and p["n"] == 180
    for c in CONDS:
        f, r = np.array(p["sig"][c]["f"]), np.array(p["sig"][c]["r"])
        k = K[c]
        assert f[50] == 40 * k and f[51] == 10 * k and f[60] == 10 * k and f[49] == 0      # sense = forward: spike at the TSS
        assert r[45:50].tolist() == [7 * k] * 5 and r[50] == 1 * k and f[45] == 0          # antisense bump upstream, in reverse
    for c in CONDS:                                    # SynM (-, 60001-70000): window starts at 55000; gene = bins 50..149
        f, r = np.array(m["sig"][c]["f"]), np.array(m["sig"][c]["r"])
        k = K[c]
        assert r[149] == 40 * k and r[148] == 10 * k and r[50] == 10 * k and r[49] == 0    # sense = REVERSE; TSS at the right end
        assert f[150:155].tolist() == [7 * k] * 5 and f[100] == 1 * k                     # antisense bump downstream (upstream of the TSS)
    assert p["strand"] == "+" and m["strand"] == "-" and p["gene_id"] == "SYNP0001.3"


def test_signal_not_rescaled_and_matches_bigwig(syn):
    data = prep(syn, loci_entries=["SynP"])
    b = np.array(data["loci"][0]["sig"]["B"]["f"])
    a = np.array(data["loci"][0]["sig"]["A"]["f"])
    assert np.allclose(b, 2 * a)                       # the bigWig values (B = 2 x A) are used as they are
    assert data["meta"]["scale_factors"] == {"A": 1.0, "B": 2.0}
    assert "already applied" in data["meta"]["normalization"]


def test_window_clipped_at_chromosome_ends_and_unknown_chromosome_reported(syn):
    data = prep(syn, loci_entries=["chrS:1-300", "chrS:599,000-600,000", "chrZ:1-100"])
    first, last = data["loci"]
    assert first["win_start"] == 1                      # 0 - flank clipped to the chromosome start
    assert last["win_start"] + last["n"] * last["bin"] - 1 <= CLEN
    assert [n["entry"] for n in data["meta"]["not_found"]] == ["chrZ:1-100"]
    assert "chrZ" in data["meta"]["not_found"][0]["reason"]


def test_wide_window_raises_bin_size_in_whole_multiples(syn):
    data = prep(syn, loci_entries=["chrS:1-590000"], flank=0, bin_size=100)
    L = data["loci"][0]
    assert L["n"] <= ss.MAX_BINS and L["bin"] % 100 == 0 and L["bin"] > 100 and L["bin_raised"]


# ---------------- gene models ----------------

def test_gene_models_from_gtf_union_of_exons_and_rows(syn):
    wins = {CHROM: [(25001, 43000)]}
    genes = ss.read_gene_models(syn["gtf"], wins)
    assert set(genes) == {"SYNP0001.3"}                                         # streamed: other genes and chrOther were never kept
    g = genes["SYNP0001.3"]
    assert g["strand"] == "+" and g["gene_name"] == "SynP" and g["extent"] == (30001, 38000)
    ext = {r.gene_id: (r.start, r.end) for r in syn["gene_table"].itertuples()}
    models, hidden = ss.models_for_window(genes, ext, CHROM, 25001, 43000)
    assert hidden == 0 and len(models) == 1
    m = models[0]
    assert m["exons"] == [[30001, 30500], [33001, 33800], [37001, 38000]]      # the extra overlapping transcript is merged away
    assert (m["start"], m["end"], m["row"]) == (30001, 38000, 0)


def test_gene_models_overlapping_genes_get_different_rows_and_cap():
    genes = {f"g{i}": {"gene_id": f"g{i}", "gene_name": f"G{i}", "chrom": "c", "strand": "+", "exons": [(100, 900)], "extent": (100, 900)}
             for i in range(13)}
    models, hidden = ss.models_for_window(genes, {}, "c", 1, 1000)
    assert len({m["row"] for m in models}) == ss.MAX_MODEL_ROWS and hidden == 3
    far = {"x": {"gene_id": "x", "gene_name": "X", "chrom": "c", "strand": "-", "exons": [(5000, 5100)], "extent": (5000, 5100)}}
    assert ss.models_for_window(far, {}, "c", 1, 1000) == ([], 0)


def test_gene_models_without_exon_lines_use_gene_extent(tmp_path):
    gtf = tmp_path / "g.gtf"
    gtf.write_text('c\tx\tgene\t100\t900\t.\t-\t.\tgene_id "g"; gene_name "G";\n')
    genes = ss.read_gene_models(str(gtf), {"c": [(1, 1000)]})
    m, _ = ss.models_for_window(genes, {}, "c", 1, 1000)
    assert m[0]["exons"] == [] and (m[0]["start"], m[0]["end"], m[0]["strand"]) == (100, 900, "-")


# ---------------- metagene ----------------

def test_orient_profiles_flips_minus_genes_and_does_not_swap_sense():
    f = np.array([[1, 2, 3], [1, 2, 3]], float)           # forward, genomic left to right
    r = np.array([[10, 20, 30], [10, 20, 30]], float)
    sint = np.array([1, -1])
    sense, anti = ss.orient_profiles(f, r, sint, True)
    assert sense[0].tolist() == [1, 2, 3] and anti[0].tolist() == [10, 20, 30]      # + gene: sense = forward, as is
    assert sense[1].tolist() == [30, 20, 10] and anti[1].tolist() == [3, 2, 1]      # - gene: sense = reverse, 5' (genomic right) first
    s2, a2 = ss.orient_profiles(f, r, sint, False)                                  # no strand column: nothing swapped or flipped
    assert s2.tolist() == f.tolist() and a2.tolist() == r.tolist()


def test_tss_metagene_5prime_left_and_antisense_not_swapped(syn):
    mg = ss.build_metagene(syn["bws"], CONDS, syn["gene_table"], ["protein_coding"], 3000, 1000, None)
    g = mg["groups"][0]
    assert g["name"] == "protein_coding" and g["n_tss"] == 25 + 1 and g["kind"] == "gene"   # 25 fillers + SynP
    v = g["tss"]
    assert len(v["x"]) == 120 and v["x"][60] == 25.0 and v["ticks"]["text"] == ["-3 kb", "TSS", "+3 kb"]
    for c in CONDS:
        k = K[c]
        p = v["profiles"][c]
        se, an = np.array(p["sense_mean"]), np.array(p["anti_mean"])
        # mean over + AND - genes: both strands line up because - genes are flipped
        assert se[:60].tolist() == [0] * 60                                    # nothing sense-side upstream of the TSS
        assert se[60] == pytest.approx(40 * k) and se[61] == pytest.approx(40 * k) and se[62] == pytest.approx(10 * k)
        assert an[50:60].tolist() == pytest.approx([7 * k] * 10) and an[60] == pytest.approx(1 * k)   # the bump is antisense
        assert not np.allclose(se, an)
        assert np.array(p["sense_median"])[60] == pytest.approx(40 * k)


def test_tss_metagene_minus_genes_only_has_spike_on_the_left(syn):
    gt = syn["gene_table"]
    minus = gt[(gt["strand"] == "-") & (gt["gene_id"].str.startswith("PC") | gt["gene_id"].str.startswith("LN"))]
    mg = ss.build_metagene(syn["bws"], ["A"], minus, ["all"], 3000, 1000, None)
    se = np.array(mg["groups"][0]["tss"]["profiles"]["A"]["sense_mean"])
    assert mg["groups"][0]["n_tss"] == len(minus) == 23
    assert se[60] == pytest.approx(40) and se[59] == 0 and se[62] == pytest.approx(10)   # 5' (TSS) on the LEFT for - genes too


def test_body_metagene_scaled_bins_and_flanks(syn):
    gt = syn["gene_table"]
    fillers = gt[gt["gene_id"].str.startswith("PC")]            # all 2 kb: 20 bp body bins (SynP's 8 kb gives 80 bp bins)
    mg = ss.build_metagene(syn["bws"], CONDS, fillers, ["protein_coding"], 3000, 1000, None)
    b = mg["groups"][0]["body"]
    assert b["ticks"]["text"] == ["-3 kb", "TSS", "TES", "+3 kb"] and len(b["x"]) == 30 + 100 + 30
    assert mg["groups"][0]["n_body"] == 25
    se = np.array(b["profiles"]["A"]["sense_mean"])
    an = np.array(b["profiles"]["A"]["anti_mean"])
    assert se[:30].tolist() == [0] * 30
    assert se[30:35].tolist() == pytest.approx([40.0] * 5)                       # first 100 bp of the 2 kb fillers = 5 body bins of 20 bp
    assert se[40] == pytest.approx(10.0)
    assert an[25:30].tolist() == pytest.approx([7.0] * 5)                          # 500 bp upstream bump = 5 flank bins of 100 bp


def test_metagene_expression_filter_and_min_group_size_and_length(syn):
    expressed = ss.expressed_gene_ids(syn["counts_m"], 10, 2)
    assert "LN20" not in expressed and "LN00" in expressed and len(expressed) == len(GENES) - 2
    assert ss.expressed_gene_ids(None, 10, 2) is None
    mg = ss.build_metagene(syn["bws"], CONDS, syn["gene_table"], ["lncRNA", "miRNA"], 3000, 1000, expressed)
    names = [g["name"] for g in mg["groups"]]
    assert names == ["lncRNA"]
    assert mg["groups"][0]["n_expressed"] == 21                                   # SynM + 22 fillers - 2 unexpressed
    assert any("'miRNA'" in n and "skipped" in n for n in mg["notes"])
    # body view needs genes >= min_length: with min_length above the filler length only SynM qualifies -> no body view
    mg2 = ss.build_metagene(syn["bws"], CONDS, syn["gene_table"], ["lncRNA"], 3000, 5000, expressed)
    assert mg2["groups"][0]["body"] is None and "fewer than 20" in mg2["groups"][0]["body_note"]


def test_lncrna_group_takes_older_gencode_types(syn):
    # GENCODE vM25 still labels many lncRNAs lincRNA / antisense; the lncRNA group must not come out empty
    expressed = ss.expressed_gene_ids(syn["counts_m"], 10, 2)
    gt = syn["gene_table"].copy()
    ln = gt.index[gt["gene_type"] == "lncRNA"]
    gt.loc[ln[::2], "gene_type"] = "lincRNA"
    gt.loc[ln[1::4], "gene_type"] = "antisense"
    mg = ss.build_metagene(syn["bws"], CONDS, gt, ["lncRNA"], 3000, 1000, expressed)
    assert [g["name"] for g in mg["groups"]] == ["lncRNA"]
    assert mg["groups"][0]["n_expressed"] == 21                                   # same genes as before the relabel
    assert "lincRNA" in mg["groups"][0]["label"] and "antisense" in mg["groups"][0]["label"]
    assert any("older GENCODE" in n for n in mg["notes"])
    # other groups still match their gene_type exactly
    only_linc = ss.build_metagene(syn["bws"], CONDS, gt, ["lincRNA"], 3000, 1000, expressed)
    assert all(g["n_expressed"] < 21 for g in only_linc["groups"])


def test_metagene_bed_oriented_and_unoriented(syn):
    # 24 centred regions: the centre of a gene body of the fillers
    regs = []
    for gid, name, typ, strand, s, e in GENES[2:27]:
        regs.append({"chrom": CHROM, "start0": s - 1, "end0": e, "name": gid, "strand": strand})
    mg = ss.build_metagene(syn["bws"], ["A"], syn["gene_table"], [], 3000, 1000, None, bed_regions=regs)
    g = mg["groups"][0]
    assert g["kind"] == "bed" and g["oriented"] and g["n_tss"] == 25 and g["body"] is None
    assert g["tss"]["ticks"]["text"][1] == "centre"
    se = np.array(g["tss"]["profiles"]["A"]["sense_mean"])
    # a 2 kb gene is centred on the middle: its 5' end (the 100 bp spike of 40) is at -1000, the body (10) runs to +1000;
    # the 1 bp reflection offset of - regions about an even-length centre makes the spike bins read ~39.6 instead of 40
    assert se[:40].tolist() == [0.0] * 40 and se[40] == pytest.approx(40, abs=1.0) and se[44:79].tolist() == pytest.approx([10.0] * 35, abs=0.5)
    assert se[80:].tolist() == pytest.approx([0.0] * 40, abs=0.2)                     # the same 1 bp offset leaks <0.2 into bin 80
    for r in regs:
        r["strand"] = ""
    mg = ss.build_metagene(syn["bws"], ["A"], syn["gene_table"], [], 3000, 1000, None, bed_regions=regs)
    g = mg["groups"][0]
    assert not g["oriented"] and any("no strand column" in n for n in mg["notes"])
    plus = np.array(g["tss"]["profiles"]["A"]["sense_mean"])                           # 'sense' = + strand RNA here, nothing swapped
    minus = np.array(g["tss"]["profiles"]["A"]["anti_mean"])
    assert plus.max() > 0 and minus.max() > 0 and not np.allclose(plus, minus)
    few = ss.build_metagene(syn["bws"], ["A"], syn["gene_table"], [], 3000, 1000, None, bed_regions=regs[:5])
    assert few["groups"] == [] and any("only 5 region(s)" in n for n in few["notes"])


def test_metagene_windows_leaving_the_chromosome_are_dropped(syn):
    gt = syn["gene_table"]
    edge = pd.DataFrame([{"gene_id": f"E{i}", "gene_name": f"E{i}", "gene_type": "protein_coding", "chr": CHROM,
                          "start": 100 + i, "end": 3100 + i, "strand": "+", "length": 3001} for i in range(25)])
    g = ss.build_metagene(syn["bws"], ["A"], edge, ["protein_coding"], 3000, 1000, None)
    assert g["groups"] == [] and any("skipped" in n or "only" in n for n in g["notes"])


def test_threads_give_the_same_profiles(syn):
    a = ss.build_metagene(syn["bws"], CONDS, syn["gene_table"], ["protein_coding"], 3000, 1000, None, threads=1)
    b = ss.build_metagene(syn["bws"], CONDS, syn["gene_table"], ["protein_coding"], 3000, 1000, None, threads=2)
    assert a["groups"][0]["tss"]["profiles"] == b["groups"][0]["tss"]["profiles"]


# ---------------- strand table ----------------

def test_strand_table_antisense_fraction_maths():
    f = pd.DataFrame({"s1": [90, 5, 0, 30], "s2": [110, 15, 0, 10], "s3": [10, 40, 0, 0]}, index=["gp", "gm", "gz", "gl"], dtype=float)
    r = pd.DataFrame({"s1": [10, 45, 0, 1], "s2": [10, 85, 0, 1], "s3": [1, 60, 0, 0]}, index=["gp", "gm", "gz", "gl"], dtype=float)
    gt = pd.DataFrame({"gene_id": ["gp", "gm", "gz", "gl"], "gene_name": ["P", "M", "Z", "L"], "gene_type": "x",
                       "chr": "c", "start": 1, "end": 2, "strand": ["+", "-", "+", "+"], "length": 2})
    cond = {"s1": "A", "s2": "A", "s3": "B"}
    df, n_pass = ss.strand_table(f, r, gt, cond, ["A", "B"], min_reads=50)
    assert n_pass == 2 and df["gene_id"].tolist() == ["gm", "gp"]                           # gz, gl below 50 reads; sorted by total reads
    gp = df.set_index("gene_id").loc["gp"]
    # + gene: sense = forward. A: sense 90+110=200, anti 10+10=20 -> 20/220 ; B: sense 10, anti 1 -> 1/11
    assert gp["A__sense"] == 200 and gp["A__antisense"] == 20
    assert gp["A__antisense_fraction"] == pytest.approx(20 / 220) and gp["B__antisense_fraction"] == pytest.approx(1 / 11)
    gm = df.set_index("gene_id").loc["gm"]
    # - gene: sense = REVERSE. A: sense 45+85=130, anti 5+15=20 -> 20/150 ; B: sense 60, anti 40 -> 0.4
    assert gm["A__sense"] == 130 and gm["A__antisense"] == 20
    assert gm["A__antisense_fraction"] == pytest.approx(20 / 150) and gm["B__antisense_fraction"] == pytest.approx(0.4)
    assert gp["total_reads"] == 200 + 20 + 10 + 1
    df0, _ = ss.strand_table(f, r, gt, cond, ["A", "B"], min_reads=0)
    assert np.isnan(df0.set_index("gene_id").loc["gz"]["A__antisense_fraction"])           # 0 / 0 -> NA, not 0
    top1, n = ss.strand_table(f, r, gt, cond, ["A", "B"], min_reads=0, top=1)
    assert len(top1) == 1 and n == 4
    block = ss.strand_table_block(df, ["A", "B"], n_pass)
    assert block["rows"][1][8][0] == [200, 20, pytest.approx(0.0909, abs=1e-4)] and block["n_shown"] == 2


def test_strand_table_in_synthetic_run(syn):
    data = prep(syn, strand_min_reads=50)
    st = data["strand_table"]
    assert st["conditions"] == CONDS and st["n_shown"] == st["n_pass"] == len(GENES) - 2        # LN20/LN21 have < 50 reads
    row = {r[0]: r for r in st["rows"]}["PC00"]       # index 2 -> sense 300 per sample, antisense 10; strand +
    assert row[6] == "+" and row[8][0] == [600, 20, pytest.approx(20 / 620, abs=1e-4)]
    rowm = {r[0]: r for r in st["rows"]}["PC01"]      # index 3, strand -: sense = reverse = 400, anti = forward = 10
    assert rowm[6] == "-" and rowm[8][1][0] == 800 and rowm[8][1][1] == 20


# ---------------- page ----------------

def strip_scripts(doc):
    return re.sub(r"<script\b.*?</script>", "", doc, flags=re.S)


def embedded_json(doc):
    m = re.search(r'<script id="ttdata" type="application/json">(.*?)</script>', doc, re.S)
    assert m
    return json.loads(m.group(1))


@pytest.fixture(scope="module")
def page(syn, tmp_path_factory):
    out = tmp_path_factory.mktemp("page") / "TTpipe_stranded_signal.html"
    bed = out.parent / "regions.bed"
    bed.write_text("chrS\t100000\t102000\tenhSyn\t0\t+\n")
    data = prep(syn, loci_entries=["SynP", "SynM", "chrS:100001-102000", "Nope"], deg_tables=syn["degs"], top_degs=1,
                regions_bed=str(bed))
    doc = ss.build_html(data)
    out.write_text(doc, encoding="utf-8")
    return {"path": out, "doc": doc, "data": data}


def test_page_self_contained_and_json_parses(page):
    doc = page["doc"]
    assert "plotly.js" in doc and not re.search(r"<script[^>]+\bsrc\s*=", doc, re.I)
    body = strip_scripts(doc)
    assert not re.search(r"<link[^>]+\bhref\s*=", body, re.I)
    assert not re.search(r"""(src|href)\s*=\s*["']\s*(https?:)?//""", body, re.I)
    ours = doc.split('<script id="ttdata"')[1].split("<script>")[1]
    assert not re.search(r"fetch\(|XMLHttpRequest|WebSocket|https?://", ours)
    emb = embedded_json(doc)
    assert [l["label"] for l in emb["loci"]][:2] == ["SynP", "SynM"]
    assert emb["meta"]["banner"] == "SYNTHETIC TEST DATA" and "SYNTHETIC TEST DATA" in body
    assert [n["entry"] for n in emb["meta"]["not_found"]] == ["Nope"] and 'id="nfbox"' in body
    assert "+ strand RNA" in body and "- strand RNA" in body
    assert "NaN" not in doc.split('id="ttdata"')[1].split("</script>")[0]


def test_embedded_data_stays_small(page):
    emb_len = len(page["doc"].split('<script id="ttdata"')[1].split("</script>")[0])
    assert emb_len < 400_000, emb_len           # binned arrays only; the plotly bundle is the rest of the file


def test_cli_end_to_end_and_static_figures(syn, tmp_path, capsys):
    out = tmp_path / "analysis" / "TTpipe_stranded_signal.html"
    args = []
    for c in CONDS:
        for s in ("forward", "reverse"):
            args += ["--bigwig", f"{c}={s}={syn['bws'][c][s]}"]
    args += ["--condition-order", "B,A", "--scale-factors", syn["sf"], "--gene-info", syn["gi"], "--gtf", syn["gtf"],
             "--counts-forward", syn["counts_f"], "--counts-reverse", syn["counts_r"], "--counts-matrix", syn["counts_m"],
             "--samples", syn["samples"], "--loci", "SynP", "SynM", "chrS:1-100", "--top-degs", "0",
             "--static-formats", "svg", "png", "--metagene-groups", "protein_coding", "--threads", "1",
             "--out", str(out)]
    assert ss.main(args) == 0
    err = capsys.readouterr().err
    assert "3 locus track(s)" in err and "protein_coding" in err and "WARNING" not in err
    emb = embedded_json(out.read_text())
    assert emb["meta"]["conditions"] == ["B", "A"]                                         # condition order respected
    files = sorted(os.listdir(tmp_path / "analysis" / "stranded_signal"))
    assert files == ["SynM.png", "SynM.svg", "SynP.png", "SynP.svg", "chrS_1-100.png", "chrS_1-100.svg"]
    svg = (tmp_path / "analysis" / "stranded_signal" / "SynP.svg").read_text()
    assert "+ strand RNA" in svg and "- strand RNA" in svg and "SynP" in svg


def test_static_formats_validated_and_matplotlib_missing_message(syn, tmp_path, monkeypatch):
    with pytest.raises(SystemExit):
        ss.main(["--bigwig", f"A=forward={syn['bws']['A']['forward']}", "--gene-info", syn["gi"], "--static-formats", "gif",
                 "--out", str(tmp_path / "x.html")])
    monkeypatch.setitem(sys.modules, "matplotlib", None)         # makes `import matplotlib` fail
    monkeypatch.setitem(sys.modules, "matplotlib.pyplot", None)
    with pytest.raises(SystemExit) as e:
        ss._import_matplotlib()
    assert "matplotlib is not installed" in str(e.value) and "track_static_formats" in str(e.value)


def test_no_choices_in_argparse():
    # regression guard for the Python 3.11 argparse trap (choices= together with nargs="*" defaults)
    src = (Path(HERE).parent / "scripts" / "stranded_signal.py").read_text()
    assert "choices=" not in src.replace("choose from", "")


def test_missing_counts_and_no_loci_still_build_a_page(syn):
    data = prep(syn, counts_forward=None, counts_reverse=None, counts_matrix=None)
    assert data["strand_table"] is None and data["loci"] == [] and data["metagene"]["groups"]
    doc = ss.build_html(data)
    assert "no strand table" in doc.lower() and 'id="noloci"' in doc


SNAKEMAKE = shutil.which("snakemake", path=os.pathsep.join([os.path.dirname(sys.executable), os.environ.get("PATH", "")]))


@pytest.mark.skipif(SNAKEMAKE is None, reason="snakemake not installed")
def test_rule_command_runs_on_synthetic_files(syn, tmp_path):
    """The real command of rule stranded_signal (workflows/rules/stranded.smk), as Snakemake renders it with `-n -p`,
    executed on the synthetic files in an analysis-only configuration with the bigWigs already present.
    (Snakemake itself is not asked to run the job: it would want conda; CONDA_PREFIX points at this environment.)"""
    root = Path(__file__).resolve().parents[1]
    out = tmp_path / "out"
    (out / "tracks").mkdir(parents=True)
    (out / "counts").mkdir()
    (out / "logs").mkdir()
    for c, new in {"A": "DMSO", "B": "2H"}.items():
        for st in ("forward", "reverse"):
            shutil.copy(syn["bws"][c][st], out / "tracks" / f"{new}_{st}.bw")
    pd.DataFrame({"condition": ["DMSO", "2H"], "mapped_pairs": [1000, 500], "scale_factor": [1.0, 2.0], "source": "computed"}).to_csv(
        out / "tracks" / "scale_factors.tsv", sep="\t", index=False)
    shutil.copy(syn["gi"], out / "counts" / "gene_info.tsv")
    ren = {"A_1": "CI_DMSO_rep1", "A_2": "CI_DMSO_rep2", "B_1": "CI_2H_rep1", "B_2": "CI_2H_rep2"}
    for src, dst in ((syn["counts_f"], "counts_forward.txt"), (syn["counts_r"], "counts_reverse.txt"), (syn["counts_m"], "counts_deseq2.txt")):
        pd.read_csv(src, sep="\t", index_col=0).rename(columns=ren).to_csv(out / "counts" / dst, sep="\t")
    target = out / "analysis" / "TTpipe_stranded_signal.html"
    cmd = [SNAKEMAKE, str(target), "-n", "-p", "-s", "workflows/TTpipe.snakefile", "--cores", "2", "--configfile",
           "config/TTpipeConfig.yaml", "tests/dryrun/config_analysis_only.yaml", "--config", f"outputDir={out}",
           f"gtf={syn['gtf']}", f"analysis_counts={out / 'counts' / 'counts_deseq2.txt'}", "track_top_degs=0",
           "track_loci=[SynP, chrS:60001-70000, Missing]", "run_analysis=False", "qc_report=True"]
    r = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-1500:] + r.stderr[-3000:]
    start = r.stdout.index("$CONDA_PREFIX/bin/python") if "$CONDA_PREFIX/bin/python" in r.stdout else r.stdout.index("/bin/python")
    lines = r.stdout[start:].splitlines()
    block = "\n".join(lines[:next(i for i, ln in enumerate(lines) if " 2> " in ln) + 1])   # the rendered command, up to the log redirect
    env = dict(os.environ, CONDA_PREFIX=sys.prefix)
    run = subprocess.run(["bash", "-c", block], cwd=root, capture_output=True, text=True, env=env)
    log_text = (out / "logs" / "stranded_signal.log").read_text()
    assert run.returncode == 0, log_text[-3000:]
    assert target.exists() and target.stat().st_size > 1_000_000
    emb = embedded_json(target.read_text())
    assert [l["label"] for l in emb["loci"]] == ["SynP", "chrS:60001-70000"]
    assert [n["entry"] for n in emb["meta"]["not_found"]] == ["Missing"]
    assert emb["meta"]["conditions"] == ["DMSO", "2H"] and emb["strand_table"]["conditions"] == ["DMSO", "2H"]
    assert {g["name"] for g in emb["metagene"]["groups"]} == {"protein_coding", "lncRNA"}
    assert "NOT FOUND: Missing" in log_text and "locus track(s) drawn" in log_text and "strand table:" in log_text
    row = {x[0]: x for x in emb["strand_table"]["rows"]}["PC00"]       # conditions come from the samplesheet names
    assert row[8][0][0] == 600 and row[8][1][0] == 600


def test_samplesheet_can_supply_the_conditions(syn, tmp_path):
    sheet = tmp_path / "sheet.txt"
    pd.DataFrame({"Name": [s for s, _ in SAMPLES], "Condition": [c for _, c in SAMPLES], "Replicate": "r"}).to_csv(
        sheet, sep="\t", index=False)
    a = prep(syn, samples=str(sheet))["strand_table"]
    b = prep(syn)["strand_table"]
    assert a == b and a["rows"][0][8][0][0] > 0
    bad = tmp_path / "bad.txt"
    bad.write_text("x\ty\n1\t2\n")
    with pytest.raises(ValueError):
        prep(syn, samples=str(bad))


def test_bigwig_spec_parsing():
    assert ss.parse_bigwigs(["A=forward=/x/a.bw", "A=reverse=/x/b.bw"]) == {"A": {"forward": "/x/a.bw", "reverse": "/x/b.bw"}}
    with pytest.raises(ValueError):
        ss.parse_bigwigs(["A=sideways=/x"])


# ---------------- the JavaScript really runs (headless Chrome) ----------------

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


def dump_dom(html_path, hash_="", extra=()):
    prof = Path(html_path).parent / "chrome_profile"
    cmd = [CHROME, "--headless=new", "--disable-gpu", "--no-sandbox", "--no-first-run", f"--user-data-dir={prof}",
           "--allow-file-access-from-files", "--virtual-time-budget=20000", *extra, "--dump-dom", f"file://{html_path}{hash_}"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        pytest.skip("headless Chrome timed out")
    if r.returncode != 0 or "<html" not in r.stdout:
        pytest.skip(f"headless Chrome unusable here (exit {r.returncode}): {r.stderr[-200:]}")
    return r.stdout


@needs_chrome
def test_javascript_runs_in_headless_chrome(page):
    f = page["path"]
    dom = dump_dom(f)
    assert 'data-ready="1"' in dom, "the page script did not finish"
    assert 'id="errbox" class="err" role="alert" hidden' in dom
    assert "js-plotly-plot" in dom and len(re.findall(r'class="point"|class="bartext"|<g class="trace bars"', dom)) > 0
    assert len(re.findall(r"<option value=\"\d+\">", dom)) >= 4            # loci in the dropdown (4: SynP, SynM, enhSyn region, top DEG)
    assert "- strand RNA" in dom
    dom = dump_dom(f, "#selftest")
    m = re.search(r'<pre id="selftest">(.*?)</pre>', dom, re.S)
    assert m, "self-test output missing"
    lines = m.group(1).strip().splitlines()
    assert lines[-1] == "SELFTEST DONE failures=0", "\n".join(lines)
    assert not [ln for ln in lines if ln.startswith("FAIL")], "\n".join(lines)
    assert sum(ln.startswith("PASS") for ln in lines) >= 25, "\n".join(lines)
    joined = "\n".join(lines)
    for needle in ("- strand traces are mirrored below zero", "y scale shared across conditions", "free y toggle",
                   "locus dropdown redraws", "search keeps matching loci only", "median toggle switches the line",
                   "metagene group switch", "gene-body view drawn", "table sorts by antisense fraction",
                   "csv rows match the shown rows"):
        assert needle in joined, needle
    print(joined)


@needs_chrome
def test_screenshot_of_the_finished_page(page, tmp_path):
    shot = Path(os.environ.get("TTPIPE_SHOT_DIR", tmp_path)) / "stranded_signal.png"
    dump_dom(page["path"], "", extra=[f"--screenshot={shot}", "--window-size=1300,1500", "--hide-scrollbars"])
    assert shot.exists() and shot.stat().st_size > 10_000
