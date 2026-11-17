"""Tests for the pure-logic R helpers and a parse check of the analysis scripts.

DESeq2 / fgsea / ggforce code paths are NOT exercised here (those packages are not required to run
the tests); only scripts/tt_analysis_utils.R (base R) is executed, plus `parse()` on every script.
"""
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not available")


def r(code, tmp_path=None):
    """Run R code with the utils sourced; return CompletedProcess."""
    prog = f'source("{SCRIPTS / "tt_analysis_utils.R"}")\n{code}\n'
    return subprocess.run(["Rscript", "--vanilla", "-e", prog], capture_output=True, text=True, cwd=tmp_path)


def r_ok(code, tmp_path=None):
    p = r(code, tmp_path)
    assert p.returncode == 0, p.stderr
    return p.stdout.strip()


@pytest.mark.parametrize("script", ["run_deseq2.R", "run_fgsea.R", "make_sankey.R", "make_rnk.R", "tt_analysis_utils.R"])
def test_scripts_parse(script):
    p = subprocess.run(["Rscript", "-e", f'invisible(parse("{SCRIPTS / script}"))'], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr


def test_parse_cli():
    out = r_ok('a <- parse_cli(c("--degs-pattern","x/{comparison}.csv","--alpha","0.05"), required="alpha", defaults=list(seed="42"));'
               'cat(a$degs_pattern, a$alpha, a$seed)')
    assert out == "x/{comparison}.csv 0.05 42"
    p = r('parse_cli(c("--a","1"), required="b")')
    assert p.returncode != 0 and "missing required" in p.stderr
    assert r('parse_cli(c("--a"))').returncode != 0


def test_validate_design():
    assert r_ok('cat(deparse(validate_design("~ condition + replicate")))') == "~condition + replicate"
    assert r("validate_design('~ condition + batch')").returncode != 0
    p = r("validate_design('~ replicate')")
    assert p.returncode != 0 and "must include 'condition'" in p.stderr
    assert r("validate_design('condition ~ replicate')").returncode != 0
    assert r("validate_design('~ condition +')").returncode != 0


def test_match_samples_by_name_not_position():
    assert r_ok('cat(match_samples(c("b","a","c"), c("a","b","c")))') == "2 1 3"
    p = r('match_samples(c("a","b"), c("a","c"))')
    assert p.returncode != 0 and "in table only: c" in p.stderr and "in counts only: b" in p.stderr
    assert r('match_samples(c("a","a"), c("a"))').returncode != 0


def test_condition_levels_and_comparisons():
    assert r_ok('cat(condition_levels(c("2H","DMSO","6H"), "DMSO"))') == "DMSO 2H 6H"
    assert r_ok('cat(condition_levels(c("2H","DMSO","6H"), "DMSO", c("DMSO","6H","2H")))') == "DMSO 6H 2H"
    assert r("condition_levels(c('a','b'), 'zzz')").returncode != 0
    assert r_ok('x <- parse_comparisons("2H:DMSO,6H:DMSO", c("DMSO","2H","6H")); cat(x[[2]], sep="|")') == "6H|DMSO"
    assert r('parse_comparisons("2H:nope", c("DMSO","2H"))').returncode != 0


def test_classify_and_summary():
    code = ('padj <- c(0.001, 0.001, 0.5, NA, 0.001, 0.01); lfc <- c(1, -1, 2, 3, NA, 1);'
            'cat(classify_deg(padj, lfc, 0.01), sep=",");'
            'cat("\\n"); s <- summarize_deg("c1", padj, lfc, 0.01); cat(s$n_up, s$n_down, s$n_tested)')
    lines = r_ok(code).splitlines()
    assert lines[0] == "Up,Down,n.s.,n.s.,n.s.,n.s."   # padj == alpha is not significant; NA padj / NA lfc -> n.s.
    assert lines[1] == "1 1 5"


def test_build_sankey_table_common_genes_and_order():
    code = '''
a <- data.frame(gene_id=c("g1","g2","g3"), log2FoldChange=c(1,-1,1), padj=c(0.001,0.001,0.9))
b <- data.frame(gene_id=c("g3","g2","g1","g4"), log2FoldChange=c(1,1,-1,1), padj=c(0.001,NA,0.001,0.001))
t <- build_sankey_table(list(c1=a, c2=b), 0.01)
write.table(t, stdout(), sep="\\t", quote=FALSE, row.names=FALSE)
t1 <- build_sankey_table(list(only=a), 0.01); cat("single:", nrow(t1), ncol(t1), "\\n")
'''
    out = r_ok(code).splitlines()
    assert out[0] == "gene_id\tc1\tc2"
    assert out[1:4] == ["g1\tUp\tDown", "g2\tDown\tn.s.", "g3\tn.s.\tUp"]   # g4 not in all comparisons
    assert out[4] == "single: 3 2"


def test_rank_by_symbol_collapse_and_sort():
    code = '''
degs <- data.frame(gene_id=c("e1","e2","e3","e4","e5","e6"), stat=c(2.0,-3.0,1.0,NA,5.0,0.5))
info <- data.frame(gene_id=c("e1","e2","e3","e4","e5","e6"), gene_name=c("A","A","B","C","","D"))
r <- rank_by_symbol(degs, info)   # A: keep e2 (|-3| > |2|); C: NA stat; e5: no symbol
write.table(r, stdout(), sep="\\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
'''
    assert r_ok(code).splitlines() == ["B\t1", "D\t0.5", "A\t-3"]


def test_rank_by_symbol_unmatched_gene_dropped():
    code = ('r <- rank_by_symbol(data.frame(gene_id=c("x","y"), stat=c(1,2)), data.frame(gene_id="y", gene_name="Y"));'
            'cat(nrow(r), r$gene_name)')
    assert r_ok(code) == "1 Y"
