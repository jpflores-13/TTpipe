#!/usr/bin/env python3
"""Scale factors for bigWig tracks: factor = (largest mapped read-pair count) / (this one's).

Input is a two-column TSV `name<TAB>mapped_pairs` (the format of the pipeline's *_final_pairs.txt
files; several files may be given). Mapped pairs = BAM records / 2, as in the guide.
Overrides (NAME=VALUE) replace the computed factor for that name (source 'config').

  scale_factors.py --pairs counts.tsv --out scale_factors.tsv [--override DMSO=1.2757 ...]
  scale_factors.py --pairs a.txt b.txt --factor-of SAMPLE      # print one factor
"""
import argparse
import sys

COLUMNS = ["condition", "mapped_pairs", "scale_factor", "source"]


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


def format_factor(x):
    return f"{x:.12g}"


def write_table(rows, handle):
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
    p.add_argument("--pairs", nargs="+", required=True)
    p.add_argument("--override", nargs="*", default=[])
    p.add_argument("--out")
    p.add_argument("--factor-of", help="print only this name's factor")
    a = p.parse_args(argv)
    rows = compute_scale_factors(read_pairs(a.pairs), parse_overrides(a.override))
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
