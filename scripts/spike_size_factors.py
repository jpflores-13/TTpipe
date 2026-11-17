#!/usr/bin/env python3
"""DESeq2 size factors / control genes from the spike-in (optional spike-in module, EXPERIMENTAL).

  --method total   size_factors.tsv (sample, size_factor) from spike_counts.tsv:
                   size_factor = spike_pairs_used / geometric_mean(spike_pairs_used over all samples).
                   Samples with more spike-in reads get a larger factor; DESeq2 divides counts by it. An explicit error is
                   raised when any sample has zero spike-in pairs (the geometric mean and the factor are undefined).
  --method genes   control_genes.txt: one prefixed spike-in gene_id per line, the spike-in genes DESeq2's
                   estimateSizeFactors(controlGenes = ...) uses. Genes are those of --spike-genes (optionally only on the
                   contigs of --use-contigs, names without prefix) that are rows of --counts-spike. At least one control
                   gene must have a positive count in every sample (median-of-ratios needs it), otherwise an error.

Standard library only; works on Python 3.8 and later.

  spike_size_factors.py --method total --spike-counts spike_counts.tsv --out size_factors.tsv
  spike_size_factors.py --method genes --spike-genes spike_genes.tsv --counts-spike counts_spike.txt --prefix spike_ \
      [--use-contigs ERCC-00002 ...] --out control_genes.txt
"""
from __future__ import annotations

import argparse
import math
import sys
from collections import OrderedDict

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import spike_summary as ss  # noqa: E402


def geometric_mean(values):
    """Geometric mean of positive numbers; ValueError for an empty list or a value <= 0."""
    values = list(values)
    if not values:
        raise ValueError("no values")
    if any(v <= 0 for v in values):
        raise ValueError("geometric mean needs values > 0")
    return math.exp(sum(math.log(v) for v in values) / len(values))


def total_size_factors(pairs):
    """OrderedDict sample -> size factor from an ordered mapping sample -> spike-in pairs."""
    zero = [s for s, n in pairs.items() if n <= 0]
    if zero:
        raise ValueError(f"sample(s) {zero} have zero spike-in read pairs (spike_pairs_used): the spike-in size factor is "
                         "undefined. Check that the spike-in was added to these samples, that the index contains the "
                         "spike-in contigs, and spike_in.use_contigs; or use spike_in.normalize: [] to only count.")
    gm = geometric_mean(pairs.values())
    return OrderedDict((s, n / gm) for s, n in pairs.items())


def read_spike_genes(path):
    """spike_genes.tsv -> OrderedDict gene_id -> contig (prefixed)."""
    out = OrderedDict()
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gi, ci = header.index("gene_id"), header.index("contig")
        for line in fh:
            if line.strip():
                f = line.rstrip("\n").split("\t")
                out[f[gi]] = f[ci]
    return out


def read_count_table(path):
    """counts_spike.txt -> (sample names, OrderedDict gene_id -> list of ints)."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = OrderedDict()
        for line in fh:
            if line.strip():
                f = line.rstrip("\n").split("\t")
                rows[f[0]] = [int(float(x)) for x in f[1:]]
    return header[1:], rows


def control_genes(gene_contigs, counts, prefix, use_contigs=()):
    """Control gene ids (in the order of `counts`): spike-in genes, on the used contigs, with a count > 0 in every sample."""
    use = set(use_contigs or ())
    genes = [g for g in counts if g in gene_contigs and (not use or gene_contigs[g][len(prefix):] in use)]
    if not genes:
        raise ValueError("no spike-in genes available as control genes: check spike_in.gtf and spike_in.use_contigs")
    ok = [g for g in genes if all(v > 0 for v in counts[g])]
    if not ok:
        raise ValueError(f"none of the {len(genes)} spike-in control gene(s) has a count > 0 in every sample, so DESeq2 "
                         "cannot estimate size factors from them. Use spike_in.size_factor_method: 'total', fewer samples "
                         "or deeper spike-in coverage.")
    return ok


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--method", required=True)
    p.add_argument("--spike-counts")
    p.add_argument("--spike-genes")
    p.add_argument("--counts-spike")
    p.add_argument("--prefix", default="spike_")
    p.add_argument("--use-contigs", nargs="*", default=[])
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)
    if a.method not in ("total", "genes"):      # checked here, not with argparse's choices option (nargs="*" defaults fail it on 3.11)
        p.error("--method must be 'total' or 'genes'")
    try:
        if a.method == "total":
            if not a.spike_counts:
                p.error("--spike-counts is required for --method total")
            rows = ss.read_counts(a.spike_counts)
            sf = total_size_factors(OrderedDict((r["sample"], r["spike_pairs_used"]) for r in rows))
            with open(a.out, "w") as fh:
                fh.write("sample\tsize_factor\n")
                for s, v in sf.items():
                    fh.write(f"{s}\t{v:.12g}\n")
        else:
            if not (a.spike_genes and a.counts_spike):
                p.error("--spike-genes and --counts-spike are required for --method genes")
            _names, counts = read_count_table(a.counts_spike)
            ok = control_genes(read_spike_genes(a.spike_genes), counts, a.prefix, a.use_contigs)
            with open(a.out, "w") as fh:
                fh.write("".join(g + "\n" for g in ok))
    except ValueError as exc:
        sys.exit(f"spike_size_factors: {exc}")


if __name__ == "__main__":
    main()
