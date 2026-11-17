#!/usr/bin/env python3
"""Merge per-sample HTSeq outputs into count matrices.

Writes counts_forward, counts_reverse (genes x samples) and the DESeq2 matrix:
  count_strand 'sum'   forward + reverse (what the guide's code does)
  count_strand 'sense' forward counts for '+' genes, reverse counts for '-' genes (gene_info strand);
                       genes whose strand is neither fall back to forward + reverse.
HTSeq summary rows (names starting '__') are dropped. Header: gene_id<TAB>sample names.

Optional spike-in mode (--spike-prefix, used only when spike_in.enabled is true; EXPERIMENTAL): the HTSeq files then
hold host genes and spike-in genes (ids starting with the prefix). The three matrices above keep the HOST genes only, in
the same order and format as without spike-in. With --spike-out the spike-in genes go to that file (gene_id x samples,
forward + reverse summed; strand-specific 'sense' counting is not applied to them because gene_info lists host genes only).
Without --spike-prefix the code path and the output are exactly those of earlier versions.
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


def split_spike(mat, prefix):
    """(host rows, spike-in rows) of a gene x sample matrix; spike-in genes are those whose id starts with `prefix`."""
    is_spike = mat.index.str.startswith(prefix)
    return mat[~is_spike], mat[is_spike]


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
    p.add_argument("--spike-prefix", default="", help="spike-in gene id prefix: spike-in genes are split off the host matrices")
    p.add_argument("--spike-out", default="", help="write the spike-in gene counts (forward + reverse) here")
    a = p.parse_args(argv)
    if a.spike_out and not a.spike_prefix:
        p.error("--spike-out needs --spike-prefix")
    fwd = build_matrix(a.forward, a.names)
    rev = build_matrix(a.reverse, a.names)
    if a.spike_prefix:
        fwd, fwd_spike = split_spike(fwd, a.spike_prefix)
        rev, rev_spike = split_spike(rev, a.spike_prefix)
        if fwd.empty:
            sys.exit(f"no host genes left after removing gene ids that start with {a.spike_prefix!r}")
        if a.spike_out:
            both = (fwd_spike + rev_spike.reindex(index=fwd_spike.index, columns=fwd_spike.columns)).astype(int)
            write_matrix(both, a.spike_out)
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
