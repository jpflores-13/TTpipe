#!/usr/bin/env python3
"""Spike-in read-pair counts per sample from `samtools idxstats` (optional spike-in module, EXPERIMENTAL).

Input: one idxstats file per sample (the final, MAPQ-filtered, deduplicated BAM), given as
`--sample NAME CONDITION PATH` (repeat for every sample). idxstats lines are
`contig<TAB>length<TAB>mapped<TAB>unmapped`; the last line `*` (unplaced unmapped reads) is ignored.

Contigs whose name starts with --prefix are spike-in contigs, all others are host contigs.
Pairs = mapped alignment records / 2, the same records/2 convention as the rest of TTpipe (scale_factors, mapq / final
pair counts); the division is applied to the summed records of each group (host, all spike-in, used spike-in), rounded
down. In the per-contig table each contig's records / 2 is rounded down on its own, so its column can sum to slightly
less than spike_pairs_all.

--use-contigs NAME ... restricts the contigs that set the normalization to these spike-in contigs (names as in the
spike-in FASTA, WITHOUT the prefix); empty = all. A name that is not a spike-in contig of the BAM is an error.

Writes
  --out-counts   spike_counts.tsv: sample, condition, host_pairs, spike_pairs_all, spike_pairs_used, spike_pct
                 spike_pct = spike_pairs_all / (host_pairs + spike_pairs_all) * 100 (blank when both are 0)
  --out-contigs  spike_contig_counts.tsv: contig (without prefix), length, used (yes/no), then one pairs column per sample

Standard library only; works on Python 3.8 and later.
"""
from __future__ import annotations

import argparse
import sys
from collections import OrderedDict

COUNT_COLUMNS = ["sample", "condition", "host_pairs", "spike_pairs_all", "spike_pairs_used", "spike_pct"]


def parse_idxstats(lines):
    """idxstats text lines -> list of (contig, length, mapped, unmapped); the `*` line and blank lines are skipped."""
    rows = []
    for n, raw in enumerate(lines, 1):
        line = raw.rstrip("\r\n")
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 4:
            raise ValueError(f"idxstats line {n}: expected 4 tab-separated columns, found {len(f)}")
        if f[0] == "*":
            continue
        try:
            rows.append((f[0], int(f[1]), int(f[2]), int(f[3])))
        except ValueError:
            raise ValueError(f"idxstats line {n}: length / mapped / unmapped are not integers: {f[1:4]}")
    return rows


def split_pairs(rows, prefix, use_contigs=()):
    """Group idxstats rows. Returns a dict with host_pairs, spike_pairs_all, spike_pairs_used and per_contig
    (OrderedDict unprefixed spike contig -> (length, pairs)). `use_contigs` are unprefixed names; empty = all."""
    use = set(use_contigs or ())
    host = spike_all = spike_used = 0
    per = OrderedDict()
    spike_names = set()
    for contig, length, mapped, _unm in rows:
        if prefix and contig.startswith(prefix):
            name = contig[len(prefix):]
            spike_names.add(name)
            per[name] = (length, mapped // 2)
            spike_all += mapped
            if not use or name in use:
                spike_used += mapped
        else:
            host += mapped
    missing = sorted(use - spike_names)
    if missing:
        raise ValueError(f"use_contigs name(s) not found among the spike-in contigs of the BAM: {missing[:10]} "
                         f"(names are those of the spike-in FASTA without the prefix {prefix!r}; "
                         f"the BAM has {len(spike_names)} spike-in contig(s))")
    return {"host_pairs": host // 2, "spike_pairs_all": spike_all // 2, "spike_pairs_used": spike_used // 2,
            "per_contig": per}


def spike_pct(host_pairs, spike_pairs_all):
    tot = host_pairs + spike_pairs_all
    return None if tot <= 0 else 100.0 * spike_pairs_all / tot


def summarise(samples, prefix, use_contigs=()):
    """samples: list of (name, condition, idxstats rows). Returns (count rows, contig table).

    count rows: list of dicts with COUNT_COLUMNS; contig table: (contig order, {contig: length}, {sample: {contig: pairs}}).
    Every sample must have the same spike-in contigs as the first."""
    out, order, lengths, per_sample = [], None, {}, OrderedDict()
    for name, cond, rows in samples:
        g = split_pairs(rows, prefix, use_contigs)
        names = list(g["per_contig"])
        if order is None:
            order = names
            lengths = {c: g["per_contig"][c][0] for c in names}
        elif set(names) != set(order):
            raise ValueError(f"sample {name!r} has a different set of spike-in contigs than {samples[0][0]!r} "
                             "(were all BAMs aligned to the same index?)")
        per_sample[name] = {c: g["per_contig"][c][1] for c in names}
        out.append({"sample": name, "condition": cond, "host_pairs": g["host_pairs"],
                    "spike_pairs_all": g["spike_pairs_all"], "spike_pairs_used": g["spike_pairs_used"],
                    "spike_pct": spike_pct(g["host_pairs"], g["spike_pairs_all"])})
    if not order:
        raise ValueError(f"no contig name starts with the spike-in prefix {prefix!r} in the idxstats files: "
                         "the BAMs were not aligned to an index that contains the spike-in "
                         "(rebuild the STAR index with spike_in.enabled: true)")
    return out, (order, lengths, per_sample)


def _fmt(v):
    return "" if v is None else (f"{v:.6g}" if isinstance(v, float) else str(v))


def write_counts(rows, handle):
    handle.write("\t".join(COUNT_COLUMNS) + "\n")
    for r in rows:
        handle.write("\t".join(_fmt(r[c]) for c in COUNT_COLUMNS) + "\n")


def write_contigs(table, prefix, use_contigs, handle):
    order, lengths, per_sample = table
    use = set(use_contigs or ())
    names = list(per_sample)
    handle.write("\t".join(["contig", "length", "used"] + names) + "\n")
    for c in order:
        used = "yes" if (not use or c in use) else "no"
        handle.write("\t".join([c, str(lengths[c]), used] + [str(per_sample[s][c]) for s in names]) + "\n")


def read_counts(path):
    """Read spike_counts.tsv -> list of dicts (numbers as int/float, blank spike_pct as None)."""
    out = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            d = dict(zip(header, line.rstrip("\n").split("\t")))
            row = {"sample": d["sample"], "condition": d["condition"]}
            for k in ("host_pairs", "spike_pairs_all", "spike_pairs_used"):
                row[k] = int(d[k])
            row["spike_pct"] = float(d["spike_pct"]) if d.get("spike_pct", "") != "" else None
            out.append(row)
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", nargs=3, action="append", metavar=("NAME", "CONDITION", "IDXSTATS"), required=True)
    p.add_argument("--prefix", default="spike_")
    p.add_argument("--use-contigs", nargs="*", default=[])
    p.add_argument("--out-counts", required=True)
    p.add_argument("--out-contigs", required=True)
    a = p.parse_args(argv)
    samples = []
    for name, cond, path in a.sample:
        with open(path) as fh:
            samples.append((name, cond, parse_idxstats(fh)))
    try:
        rows, table = summarise(samples, a.prefix, a.use_contigs)
    except ValueError as exc:
        sys.exit(f"spike_summary: {exc}")
    with open(a.out_counts, "w") as fh:
        write_counts(rows, fh)
    with open(a.out_contigs, "w") as fh:
        write_contigs(table, a.prefix, a.use_contigs, fh)
    for r in rows:
        if r["spike_pairs_used"] == 0:
            print(f"spike_summary: warning: sample {r['sample']} has no spike-in read pairs in the contigs used",
                  file=sys.stderr)
    print(f"spike_summary: {len(rows)} sample(s), {len(table[0])} spike-in contig(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
