#!/usr/bin/env python3
"""Scale factors for bigWig tracks: factor = (largest mapped read-pair count) / (this one's).

Input is a two-column TSV `name<TAB>mapped_pairs` (the format of the pipeline's *_final_pairs.txt
files; several files may be given). Mapped pairs = BAM records / 2, as in the guide.
Overrides (NAME=VALUE) replace the computed factor for that name (source 'config').

  scale_factors.py --pairs counts.tsv --out scale_factors.tsv [--override DMSO=1.2757 ...]
  scale_factors.py --pairs a.txt b.txt --factor-of SAMPLE      # print one factor

Optional spike-in mode (EXPERIMENTAL; only used when spike_in.enabled is true in the config): the pair counts come from
spike_counts.tsv (scripts/spike_summary.py) instead of --pairs, summed over the replicates of each condition
(--group-by condition) or per sample (--group-by sample):
  --basis spike   factor = (largest group's spike_pairs_used) / (this group's), source 'spike-in'
  --basis host    factor = (largest group's host_pairs) / (this group's), source 'computed'; host_pairs exclude the
                  spike-in contigs, so the factor is the library-size factor of the host reads only
The table then has the extra columns host_pairs and spike_pairs (spike_pairs_used); mapped_pairs is the count the factor is
computed from. An override still wins (source 'config'). A group with zero spike-in pairs is an error.
  scale_factors.py --spike-counts spike_counts.tsv --group-by condition --basis spike --order DMSO 2H --out scale_factors.tsv
"""
import argparse
import os
import sys

COLUMNS = ["condition", "mapped_pairs", "scale_factor", "source"]
SPIKE_COLUMNS = COLUMNS + ["host_pairs", "spike_pairs"]


def read_pairs(paths):
    """Ordered dict name -> pairs from TSV files; a repeated name is an error."""
    pairs = {}
    for path in paths:
        with open(path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                name, n = line.rstrip("\n").split("\t")[:2]
                if name in pairs:
                    raise ValueError(f"{name!r} appears twice in the pair counts")
                pairs[name] = int(n)
    return pairs


def compute_scale_factors(pairs, overrides=None):
    """Return rows (name, pairs, factor, source) in the order of `pairs`.

    factor = max(pairs) / pairs[name]; an entry in `overrides` wins (source 'config').
    The maximum is taken over ALL names, including overridden ones.
    """
    overrides = {str(k): float(v) for k, v in (overrides or {}).items()}
    unknown = set(overrides) - set(pairs)
    if unknown:
        raise ValueError(f"scale_factors override(s) for unknown condition(s): {sorted(unknown)}")
    if not pairs:
        raise ValueError("no pair counts given")
    if any(n <= 0 for n in pairs.values()):
        bad = sorted(k for k, n in pairs.items() if n <= 0)
        if not all(k in overrides for k in bad):
            raise ValueError(f"zero mapped pairs for {bad}: cannot compute a scale factor")
    top = max(pairs.values())
    rows = []
    for name, n in pairs.items():
        if name in overrides:
            rows.append((name, n, overrides[name], "config"))
        else:
            rows.append((name, n, top / n, "computed"))
    return rows


def group_spike_counts(rows, group_by="condition", order=None):
    """Sum host_pairs / spike_pairs_used of spike_counts rows (dicts) per condition or per sample.

    Returns an ordered dict name -> (host_pairs, spike_pairs); `order` (condition names) fixes the order."""
    if group_by not in ("condition", "sample"):
        raise ValueError(f"group_by must be 'condition' or 'sample', got {group_by!r}")
    groups = {}
    for r in rows:
        g = groups.setdefault(r[group_by], [0, 0])
        g[0] += int(r["host_pairs"])
        g[1] += int(r["spike_pairs_used"])
    names = list(order) if order else list(groups)
    unknown = [n for n in names if n not in groups]
    if unknown:
        raise ValueError(f"{group_by}(s) {unknown} are not in the spike-in counts")
    return {n: tuple(groups[n]) for n in names}


def compute_spike_scale_factors(groups, basis="spike", overrides=None):
    """Rows (name, mapped_pairs, factor, source, host_pairs, spike_pairs) from `groups` (name -> (host, spike) pairs).

    basis 'spike': factor = max(spike) / spike[name], source 'spike-in'; basis 'host': factor = max(host) / host[name],
    source 'computed'. An override wins (source 'config'); the maximum is taken over all groups, overridden ones included.
    """
    if basis not in ("spike", "host"):
        raise ValueError(f"basis must be 'spike' or 'host', got {basis!r}")
    overrides = {str(k): float(v) for k, v in (overrides or {}).items()}
    unknown = set(overrides) - set(groups)
    if unknown:
        raise ValueError(f"scale_factors override(s) for unknown condition(s): {sorted(unknown)}")
    if not groups:
        raise ValueError("no pair counts given")
    idx = 1 if basis == "spike" else 0
    counts = {k: v[idx] for k, v in groups.items()}
    bad = sorted(k for k, n in counts.items() if n <= 0 and k not in overrides)
    if bad:
        if basis == "spike":
            raise ValueError(f"zero spike-in read pairs for {bad}: cannot compute a spike-in scale factor. Check that the "
                             "spike-in was added to these samples, that the STAR index contains the spike-in contigs, "
                             "and spike_in.use_contigs; or remove 'tracks' from spike_in.normalize.")
        raise ValueError(f"zero mapped host pairs for {bad}: cannot compute a scale factor")
    top = max(counts.values())
    rows = []
    for name, n in counts.items():
        host, spike = groups[name]
        if name in overrides:
            rows.append((name, n, overrides[name], "config", host, spike))
        else:
            rows.append((name, n, top / n, "spike-in" if basis == "spike" else "computed", host, spike))
    return rows


def format_factor(x):
    return f"{x:.12g}"


def write_table(rows, handle):
    if rows and len(rows[0]) > 4:                  # spike-in mode: host_pairs, spike_pairs
        handle.write("\t".join(SPIKE_COLUMNS) + "\n")
        for name, n, sf, src, host, spike in rows:
            handle.write(f"{name}\t{n}\t{format_factor(sf)}\t{src}\t{host}\t{spike}\n")
        return
    handle.write("\t".join(COLUMNS) + "\n")
    for name, n, sf, src in rows:
        handle.write(f"{name}\t{n}\t{format_factor(sf)}\t{src}\n")


def parse_overrides(items):
    out = {}
    for it in items or []:
        k, sep, v = it.rpartition("=")
        if not sep or not k:
            raise ValueError(f"override must be NAME=VALUE, got {it!r}")
        out[k] = float(v)
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pairs", nargs="+")
    p.add_argument("--override", nargs="*", default=[])
    p.add_argument("--out")
    p.add_argument("--factor-of", help="print only this name's factor")
    p.add_argument("--spike-counts", help="spike_counts.tsv: spike-in mode (see the module docstring)")
    p.add_argument("--group-by", default="condition")
    p.add_argument("--basis", default="spike")
    p.add_argument("--order", nargs="*", default=[], help="spike-in mode: condition order of the table")
    a = p.parse_args(argv)
    if a.spike_counts:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import spike_summary
        try:
            groups = group_spike_counts(spike_summary.read_counts(a.spike_counts), a.group_by, a.order or None)
            rows = compute_spike_scale_factors(groups, a.basis, parse_overrides(a.override))
        except ValueError as exc:
            sys.exit(f"scale_factors: {exc}")
    elif a.pairs:
        rows = compute_scale_factors(read_pairs(a.pairs), parse_overrides(a.override))
    else:
        p.error("one of --pairs or --spike-counts is required")
    if a.factor_of:
        hit = [r for r in rows if r[0] == a.factor_of]
        if not hit:
            sys.exit(f"{a.factor_of!r} not in the pair counts")
        print(format_factor(hit[0][2]))
    if a.out:
        with open(a.out, "w") as fh:
            write_table(rows, fh)
    elif not a.factor_of:
        write_table(rows, sys.stdout)


if __name__ == "__main__":
    main()
