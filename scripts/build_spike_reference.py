#!/usr/bin/env python3
"""Combined host + spike-in reference for the optional spike-in module (EXPERIMENTAL).

Writes, from the host genome FASTA / GTF and the spike-in FASTA / (optional) GTF:

  --out-fasta   host FASTA followed by the spike-in FASTA, every spike-in header prefixed (e.g. ERCC-00002 -> spike_ERCC-00002)
  --out-gtf     host GTF followed by the spike-in GTF with seqname, gene_id and transcript_id prefixed
  --out-contigs spike_contigs.tsv: header `contig<TAB>length`, one row per prefixed spike-in contig
  --out-genes   (only with --spike-gtf) spike_genes.tsv: gene_id, contig, start, end, strand of every spike-in gene

Aligning to ONE combined reference means a read can only align to the host or to the spike-in, never to both.

Spike-in `gene` lines. HTSeq counts `--type gene`. A spike-in GTF that has only `exon` (and/or `transcript`) lines gets one
synthesised `gene` line per gene_id that has none (start = smallest start, end = largest end, strand and seqname of the
first line seen for that gene_id; source `TTpipe`, attributes `gene_id` and `gene_type "spike_in"`). Genes that already
have a `gene` line are left as they are.

Checks (the run stops with a clear message): a host contig or host GTF seqname that already starts with the prefix
(spike-in and host reads could no longer be told apart by name); a spike-in contig name used twice; a spike-in GTF
seqname that is not a contig of the spike-in FASTA; a header without a name; a contig without sequence; a spike-in
GTF line that has a gene-level feature but no gene_id.

All inputs may be plain or gzip-compressed (.gz). Files are read and written line by line, so memory use does not
depend on the genome size. Standard library only; works on Python 3.8 and later.

  build_spike_reference.py --host-fasta genome.fa --host-gtf genes.gtf --spike-fasta ERCC92.fa [--spike-gtf ERCC92.gtf]
      --prefix spike_ --out-fasta host_spike.fa --out-gtf host_spike.gtf --out-contigs spike_contigs.tsv [--out-genes spike_genes.tsv]
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import OrderedDict

PREFIX_OK = re.compile(r"[A-Za-z0-9._\-]+")
# one GTF attribute: key, whitespace, "quoted value" or bare value, optional ';' (same grammar as gene_info.py)
_ATTR = re.compile(r'(\s*)([^\s";]+)(\s+)(?:"([^"]*)"|([^;"]*?))(\s*;?)(?=\s|$|;)')
_ID_KEYS = ("gene_id", "transcript_id")


class ReferenceError(Exception):
    """A problem with the input files that the user can fix; the message says what and where."""


def open_text(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def _clean(line):
    return line.rstrip("\r\n")


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------

def copy_fasta(lines, out, prefix="", label="FASTA"):
    """Copy FASTA `lines` to `out` (blank lines dropped, '\\n' endings), prepending `prefix` to every header name.

    Returns an OrderedDict contig name (as written, i.e. with the prefix) -> sequence length. Raises ReferenceError for
    sequence before the first header, an empty header name, a repeated name, or a contig without sequence.
    """
    contigs = OrderedDict()
    current = None
    for n, raw in enumerate(lines, 1):
        line = _clean(raw)
        if not line.strip():
            continue
        if line.startswith(">"):
            body = line[1:]
            parts = body.split(None, 1)
            if not parts:
                raise ReferenceError(f"{label} line {n}: a '>' header without a sequence name")
            name = prefix + parts[0]
            if name in contigs:
                raise ReferenceError(f"{label}: the sequence name {parts[0]!r} (line {n}) is used twice")
            contigs[name] = 0
            current = name
            out.write(">" + name + ((" " + parts[1]) if len(parts) > 1 else "") + "\n")
        else:
            if current is None:
                raise ReferenceError(f"{label} line {n}: sequence before the first '>' header (not a FASTA file?)")
            seq = line.strip()
            contigs[current] += len(seq)
            out.write(seq + "\n")
    empty = [c for c, n in contigs.items() if n == 0]
    if empty:
        raise ReferenceError(f"{label}: sequence(s) without any bases: {', '.join(empty[:5])}")
    return contigs


# --------------------------------------------------------------------------
# GTF
# --------------------------------------------------------------------------

def parse_attributes(text):
    """GTF attribute column -> dict (first value of a key wins)."""
    attrs = {}
    for m in _ATTR.finditer(text):
        val = m.group(4) if m.group(4) is not None else (m.group(5) or "").strip()
        attrs.setdefault(m.group(2), val)
    return attrs


def prefix_attributes(text, prefix):
    """Prefix the values of gene_id and transcript_id in a GTF attribute column; every other character is kept.

    The prefixed values are always written in double quotes."""
    def repl(m):
        if m.group(2) not in _ID_KEYS:
            return m.group(0)
        val = m.group(4) if m.group(4) is not None else (m.group(5) or "").strip()
        return f'{m.group(1)}{m.group(2)}{m.group(3)}"{prefix}{val}"{m.group(6)}'
    return _ATTR.sub(repl, text)


def process_spike_gtf(lines, out, prefix, contigs, label="spike-in GTF"):
    """Write the prefixed spike-in GTF lines to `out`; return (genes, synthesised).

    `contigs` are the prefixed contig names of the spike-in FASTA. `genes` is an OrderedDict
    prefixed gene_id -> {contig, start, end, strand, has_gene_line}; `synthesised` lists the gene ids that got a
    synthesised `gene` line (the lines are appended to `out` after all other lines).
    """
    genes = OrderedDict()
    bad_seq = []
    for n, raw in enumerate(lines, 1):
        line = _clean(raw)
        if not line.strip():
            continue
        if line.startswith("#"):
            out.write(line + "\n")
            continue
        f = line.split("\t")
        if len(f) < 9:
            raise ReferenceError(f"{label} line {n}: expected 9 tab-separated columns, found {len(f)}")
        seq = prefix + f[0]
        if seq not in contigs:
            bad_seq.append(f[0])
            continue
        try:
            start, end = int(f[3]), int(f[4])
        except ValueError:
            raise ReferenceError(f"{label} line {n}: start/end are not integers ({f[3]!r}, {f[4]!r})")
        attrs = parse_attributes(f[8])
        gid = attrs.get("gene_id")
        if not gid:
            raise ReferenceError(f"{label} line {n} ({f[2]}) has no gene_id attribute; every spike-in GTF line needs one")
        pgid = prefix + gid
        g = genes.get(pgid)
        if g is None:
            g = genes[pgid] = {"contig": seq, "start": start, "end": end, "strand": f[6], "has_gene_line": False}
        else:
            g["start"], g["end"] = min(g["start"], start), max(g["end"], end)
        if f[2] == "gene":
            g["has_gene_line"] = True
            g["start"], g["end"], g["strand"], g["contig"] = start, end, f[6], seq
        f[0] = seq
        f[8] = prefix_attributes(f[8], prefix)
        out.write("\t".join(f) + "\n")
    if bad_seq:
        raise ReferenceError(f"{label}: sequence name(s) {sorted(set(bad_seq))[:5]} are not in the spike-in FASTA "
                             "(names are compared without the prefix)")
    synth = []
    for pgid, g in genes.items():
        if g["has_gene_line"]:
            continue
        synth.append(pgid)
        out.write("\t".join([g["contig"], "TTpipe", "gene", str(g["start"]), str(g["end"]), ".", g["strand"], ".",
                             f'gene_id "{pgid}"; gene_type "spike_in";']) + "\n")
    return genes, synth


def host_gtf_check(lines, out, prefix):
    """Copy the host GTF; return the set of its seqnames that start with the prefix (should be empty)."""
    clash = set()
    for raw in lines:
        line = _clean(raw)
        if not line.strip():
            continue
        out.write(line + "\n")
        if not line.startswith("#"):
            seq = line.split("\t", 1)[0]
            if prefix and seq.startswith(prefix):
                clash.add(seq)
    return clash


# --------------------------------------------------------------------------
# Whole build
# --------------------------------------------------------------------------

def build(host_fasta, host_gtf, spike_fasta, spike_gtf, prefix, out_fasta, out_gtf, out_contigs, out_genes=None):
    """Write the combined reference. Returns a summary dict (counts only); raises ReferenceError on bad input."""
    if not prefix or not PREFIX_OK.fullmatch(prefix):
        raise ReferenceError(f"prefix {prefix!r} must be non-empty and contain only letters, digits, '.', '_' and '-'")
    if out_genes and not spike_gtf:
        raise ReferenceError("--out-genes needs --spike-gtf")
    # spike-in FASTA first into memory-free counters: we need its contig set before the GTF, and its lengths at the end
    with open_text(host_fasta) as hf, open(out_fasta, "w") as fa:
        host_contigs = copy_fasta(hf, fa, "", "host FASTA")
        clash = sorted(c for c in host_contigs if c.startswith(prefix))
        if clash:
            raise ReferenceError(
                f"the host genome has contig(s) whose name starts with the spike-in prefix {prefix!r}: {', '.join(clash[:5])}"
                " (spike-in and host reads could no longer be told apart). Choose another spike_in.prefix.")
        with open_text(spike_fasta) as sf:
            spike_contigs = copy_fasta(sf, fa, prefix, "spike-in FASTA")
    if not spike_contigs:
        raise ReferenceError(f"no sequences found in the spike-in FASTA {spike_fasta}")
    with open_text(host_gtf) as hg, open(out_gtf, "w") as gtf:
        gclash = host_gtf_check(hg, gtf, prefix)
        if gclash:
            raise ReferenceError(f"the host GTF has seqname(s) starting with the spike-in prefix {prefix!r}: "
                                 f"{', '.join(sorted(gclash)[:5])}. Choose another spike_in.prefix.")
        genes, synth = OrderedDict(), []
        if spike_gtf:
            with open_text(spike_gtf) as sg:
                genes, synth = process_spike_gtf(sg, gtf, prefix, spike_contigs)
    with open(out_contigs, "w") as fh:
        fh.write("contig\tlength\n")
        for c, n in spike_contigs.items():
            fh.write(f"{c}\t{n}\n")
    if out_genes:
        with open(out_genes, "w") as fh:
            fh.write("gene_id\tcontig\tstart\tend\tstrand\n")
            for gid, g in genes.items():
                fh.write(f"{gid}\t{g['contig']}\t{g['start']}\t{g['end']}\t{g['strand']}\n")
    return {"host_contigs": len(host_contigs), "spike_contigs": len(spike_contigs), "spike_genes": len(genes),
            "synthesised_genes": len(synth)}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--host-fasta", required=True)
    p.add_argument("--host-gtf", required=True)
    p.add_argument("--spike-fasta", required=True)
    p.add_argument("--spike-gtf", default="")
    p.add_argument("--prefix", default="spike_")
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-gtf", required=True)
    p.add_argument("--out-contigs", required=True)
    p.add_argument("--out-genes", default="")
    a = p.parse_args(argv)
    try:
        s = build(a.host_fasta, a.host_gtf, a.spike_fasta, a.spike_gtf or None, a.prefix, a.out_fasta, a.out_gtf,
                  a.out_contigs, a.out_genes or None)
    except (ReferenceError, OSError) as exc:
        sys.exit(f"build_spike_reference: {exc}")
    print(f"build_spike_reference: {s['host_contigs']} host contig(s) + {s['spike_contigs']} spike-in contig(s); "
          f"{s['spike_genes']} spike-in gene(s), {s['synthesised_genes']} with a synthesised gene line", file=sys.stderr)


if __name__ == "__main__":
    main()
