#!/usr/bin/env python3
"""Merge per-sample HTSeq outputs into count matrices.

Writes counts_forward, counts_reverse (genes x samples) and the DESeq2 matrix:
  count_strand 'sum'   forward + reverse (what the guide's code does)
  count_strand 'sense' forward counts for '+' genes, reverse counts for '-' genes (gene_info strand);
                       genes whose strand is neither fall back to forward + reverse.
HTSeq summary rows (names starting '__') are dropped. Header: gene_id<TAB>sample names.
"""
import argparse
import sys

import pandas as pd


def read_htseq(path):
    """One HTSeq file -> int Series indexed by gene id, `__` summary rows removed."""
    df = pd.read_csv(path, sep="\t", header=None, names=["gene_id", "count"], dtype={"gene_id": str},
                     keep_default_na=False)
    df = df[~df["gene_id"].str.startswith("__")]
    if df["gene_id"].duplicated().any():
        raise ValueError(f"{path}: duplicated gene ids")
    return df.set_index("gene_id")["count"].astype(int)


def build_matrix(paths, names):
    """Columns in the order of `names`; rows in the order of the first file. Gene sets must agree."""
    if len(paths) != len(names):
        raise ValueError("number of HTSeq files and sample names differ")
    cols = {}
    for p, n in zip(paths, names):
        cols[n] = read_htseq(p)
    first = cols[names[0]].index
    for n in names[1:]:
        if set(cols[n].index) != set(first):
            raise ValueError(f"sample {n!r} has a different gene set than {names[0]!r}")
    mat = pd.DataFrame({n: cols[n].reindex(first) for n in names})
    mat.index.name = "gene_id"
    return mat.astype(int)


def deseq2_matrix(fwd, rev, mode, gene_info=None):
    if not fwd.index.equals(rev.index) or list(fwd.columns) != list(rev.columns):
        rev = rev.reindex(index=fwd.index, columns=fwd.columns)
    if mode == "sum":
        return (fwd + rev).astype(int)
    if mode != "sense":
        raise ValueError(f"count_strand must be 'sum' or 'sense', got {mode!r}")
    strand = gene_info["strand"].reindex(fwd.index)
    if strand.isna().any():
        missing = list(strand.index[strand.isna()][:5])
        raise ValueError(f"{int(strand.isna().sum())} counted genes are missing from gene_info, e.g. {missing}")
    import numpy as np
    plus = np.broadcast_to((strand == "+").to_numpy()[:, None], fwd.shape)
    minus = np.broadcast_to((strand == "-").to_numpy()[:, None], fwd.shape)
    vals = np.where(plus, fwd.to_numpy(), np.where(minus, rev.to_numpy(), fwd.to_numpy() + rev.to_numpy()))
    out = pd.DataFrame(vals, index=fwd.index, columns=fwd.columns)
    return out.astype(int)


def read_gene_info(path):
    return pd.read_csv(path, sep="\t", dtype={"gene_id": str, "strand": str}, keep_default_na=False).set_index("gene_id")


def write_matrix(df, path):
    df.to_csv(path, sep="\t", index=True, index_label="gene_id")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--names", nargs="+", required=True, help="sample names, in output column order")
    p.add_argument("--forward", nargs="+", required=True, help="HTSeq files, same order as --names")
    p.add_argument("--reverse", nargs="+", required=True)
    p.add_argument("--count-strand", choices=["sum", "sense"], default="sum")
    p.add_argument("--gene-info")
    p.add_argument("--out-forward", required=True)
    p.add_argument("--out-reverse", required=True)
    p.add_argument("--out-deseq2", required=True)
    a = p.parse_args(argv)
    fwd = build_matrix(a.forward, a.names)
    rev = build_matrix(a.reverse, a.names)
    gi = None
    if a.count_strand == "sense":
        if not a.gene_info:
            sys.exit("--gene-info is required for count_strand 'sense'")
        gi = read_gene_info(a.gene_info)
    des = deseq2_matrix(fwd, rev, a.count_strand, gi)
    write_matrix(fwd, a.out_forward)
    write_matrix(rev, a.out_reverse)
    write_matrix(des, a.out_deseq2)


if __name__ == "__main__":
    main()
