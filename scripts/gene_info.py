#!/usr/bin/env python3
"""Gene table from the `gene` lines of a GTF (plain or .gz).

Columns: gene_id, gene_name, gene_type, chr, start, end, strand, length
  - gene_type comes from `gene_type` (GENCODE) or `gene_biotype` (Ensembl), else empty.
  - gene_name falls back to the id when the GTF has no gene_name.
  - start/end are the GTF's 1-based inclusive coordinates; length = end - start + 1.
  - a repeated id keeps its first line.
"""
import argparse
import gzip
import re
import sys

COLUMNS = ["gene_id", "gene_name", "gene_type", "chr", "start", "end", "strand", "length"]
_ATTR = re.compile(r'\s*([^\s"]+)\s+(?:"([^"]*)"|([^;"]*))\s*;?')


def parse_attributes(text):
    """Parse a GTF attribute column (key "value"; key "value";) into a dict (first value wins)."""
    attrs = {}
    for m in _ATTR.finditer(text):
        key = m.group(1)
        val = m.group(2) if m.group(2) is not None else (m.group(3) or "").strip()
        attrs.setdefault(key, val)
    return attrs


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_gene_info(lines, id_attribute="gene_id"):
    """Return a list of row dicts (see COLUMNS) from an iterable of GTF lines."""
    rows, seen = [], set()
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        a = parse_attributes(f[8])
        gid = a.get(id_attribute)
        if not gid or gid in seen:
            continue
        seen.add(gid)
        start, end = int(f[3]), int(f[4])
        rows.append({
            "gene_id": gid,
            "gene_name": a.get("gene_name") or gid,
            "gene_type": a.get("gene_type") or a.get("gene_biotype") or "",
            "chr": f[0], "start": start, "end": end, "strand": f[6],
            "length": end - start + 1,
        })
    return rows


def write_gene_info(rows, handle):
    handle.write("\t".join(COLUMNS) + "\n")
    for r in rows:
        handle.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gtf", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--id-attribute", default="gene_id", help="GTF attribute used as gene_id (default gene_id)")
    a = p.parse_args(argv)
    with _open(a.gtf) as fh:
        rows = read_gene_info(fh, a.id_attribute)
    if not rows:
        sys.exit(f"no 'gene' lines with attribute {a.id_attribute!r} found in {a.gtf}")
    with open(a.out, "w") as out:
        write_gene_info(rows, out)


if __name__ == "__main__":
    main()
