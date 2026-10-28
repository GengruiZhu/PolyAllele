#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qc_matrix.py -- quality control for a PolyAllele allele matrix.

    python bin/qc_matrix.py --matrix polyallele_out/polyallele.allele_matrix.tsv
    python bin/qc_matrix.py --matrix M.tsv --bed-dir work/target --report qc.txt

Reports, per matrix:

  coverage            genes placed vs genes in the annotation (if a GFF is given)
  subgenome occupancy how many columns each row fills  -> should peak near the
                      chromosome count of the homology group
  1:1 rate            cells holding exactly one gene
  copy-number spread  genes per cell; extras are real duplications, not errors
  span purity         for multi-gene cells, the physical distance between the
                      genes -- requires --bed-dir
  over-split check    rows sharing the same Ref_Gene (should be 1)

Exit code is always 0; this is a report, not a gate.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
from collections import Counter, defaultdict

import numpy as np

NULLS = {"NA", "NAN", "", ".", "-"}


def split_cell(v):
    v = str(v).strip()
    if v.upper() in NULLS:
        return []
    return [g.strip() for g in v.split(",") if g.strip() and g.strip().upper() not in NULLS]


def load_beds(bed_dir):
    pos = {}
    if not bed_dir:
        return pos
    for bed in glob.glob(os.path.join(bed_dir, "**", "*.bed"), recursive=True):
        with open(bed) as fh:
            for i, line in enumerate(fh):
                p = line.rstrip("\n").split("\t")
                if len(p) >= 4:
                    try:
                        pos[p[3]] = (p[0], int(p[1]), int(p[2]), i)
                    except ValueError:
                        continue
    return pos


def count_genes_gff(path):
    n = 0
    if not path or not os.path.exists(path):
        return 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split("\t")
            if len(p) >= 3 and p[2] == "gene":
                n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description="QC report for a PolyAllele matrix")
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--bed-dir", default="", help="target BED dir, for span purity")
    ap.add_argument("--gff", default="", help="target GFF3, for coverage")
    ap.add_argument("--report", default="", help="also write the text to this file")
    args = ap.parse_args()

    out = []

    def say(*a):
        line = " ".join(str(x) for x in a)
        print(line)
        out.append(line)

    pos = load_beds(args.bed_dir)
    n_gff = count_genes_gff(args.gff)

    cells = []                     # genes per (row, column)
    row_sizes = []
    occ = Counter()
    synthetic_rows = 0
    ref2rows = defaultdict(set)
    gene_seen = Counter()
    n_rows = 0
    multi_cells_genes = []         # (list of genes) for multi-gene cells

    with open(args.matrix, newline="") as fh:
        r = csv.reader(fh, delimiter="\t")
        head = next(r)
        ai = [(c, i) for i, c in enumerate(head) if c.startswith("Allele_")]
        ei = head.index("Extra_Copies") if "Extra_Copies" in head else None
        gi = head.index("HomologyGroup") if "HomologyGroup" in head else None
        ri = head.index("Ref_Gene") if "Ref_Gene" in head else None
        si = head.index("Synthetic") if "Synthetic" in head else None
        for row in r:
            n_rows += 1
            syn = (si is not None and row[si] == "1")
            synthetic_rows += 1 if syn else 0
            per = defaultdict(list)
            for name, i in ai:
                if i < len(row):
                    g = split_cell(row[i])
                    if g:
                        per[name[-1]].extend(g)
            if ei is not None and ei < len(row):
                for tok in split_cell(row[ei]):
                    if ":" in tok:
                        c, g = tok.rsplit(":", 1)
                        per[c[-1]].append(g)
                    else:
                        per["?"].append(tok)
            occ[len(per)] += 1
            tot = 0
            for col, gs in per.items():
                cells.append(len(gs))
                tot += len(gs)
                for g in gs:
                    gene_seen[g] += 1
                if len(gs) > 1:
                    multi_cells_genes.append(gs)
            row_sizes.append(tot)
            if ri is not None and ri < len(row) and row[ri] not in NULLS and not syn:
                ref2rows[row[ri]].add(row[0])

    cells = np.array(cells) if cells else np.array([0])
    rs = np.array([x for x in row_sizes if x > 0]) if any(row_sizes) else np.array([0])
    placed = int(cells.sum())

    say("=" * 88)
    say(" PolyAllele -- matrix QC report")
    say(f" matrix: {os.path.abspath(args.matrix)}")
    say("=" * 88)
    say("")
    say("[ 1 ] SIZE")
    say(f"  rows (loci)                    : {n_rows:,}")
    say(f"    of which synthetic-interval  : {synthetic_rows:,} "
        f"({100.0*synthetic_rows/max(n_rows,1):.1f}%)")
    say(f"  filled cells                   : {len(cells):,}")
    say(f"  genes placed                   : {placed:,}")
    if n_gff:
        say(f"  genes in the annotation        : {n_gff:,}")
        say(f"  coverage                       : {100.0*placed/n_gff:.2f}%")
    say("")
    say("[ 2 ] COPIES PER CELL   (duplication WITHIN one subgenome)")
    say("  a cell is (locus, chromosome).  More than one gene in it means that")
    say("  subgenome carries several copies at that locus -- tandem, proximal or")
    say("  dispersed.  This is NOT the same as 'the gene is single copy':")
    say("    * a gene with no paralogue anywhere is single copy in the usual sense")
    say("    * a row with only one filled cell is a COVERAGE statement -- the other")
    say("      subgenomes' copies may be absent, or may simply not have been placed")
    say("    * 'one gene in every column' is a complete 1:1 homology group")
    n1 = int((cells == 1).sum())
    say(f"  cells with exactly 1 gene      : {n1:,} "
        f"({100.0*n1/max(len(cells),1):.2f}%)")
    say(f"  multi-gene cells               : {len(cells)-n1:,} "
        f"({100.0*(len(cells)-n1)/max(len(cells),1):.2f}%)")
    say(f"  inflation (genes / cells)      : {placed/max(len(cells),1):.4f}x")
    say("  genes per cell tail:")
    cc = Counter(cells.tolist())
    for k in sorted(cc)[:8]:
        say(f"      {k:3d} gene(s): {cc[k]:9,}  ({100.0*cc[k]/len(cells):6.3f}%)")
    for thr in (2, 5, 10):
        say(f"      cells with > {thr:2d} gene(s): {sum(v for k,v in cc.items() if k>thr):,}")
    say("")
    say("[ 3 ] SUBGENOME OCCUPANCY   (columns filled per row)")
    say("  a healthy polyploid matrix peaks near the chromosome count of the group")
    tot = sum(occ.values())
    for k in sorted(occ):
        say(f"      {k:3d} col(s): {occ[k]:9,}  ({100.0*occ[k]/max(tot,1):6.2f}%)")
    say("")
    say("[ 4 ] ROW SIZE")
    say(f"  mean={rs.mean():.2f}  median={np.median(rs):.0f}  "
        f"p90={np.percentile(rs,90):.0f}  p99={np.percentile(rs,99):.0f}  "
        f"max={rs.max():.0f}")
    say(f"  singleton rows                 : {int((rs==1).sum()):,} "
        f"({100.0*(rs==1).sum()/max(len(rs),1):.1f}%)")
    say("")
    say("[ 5 ] OVER-MERGE CHECK   (one reference locus must occupy exactly one row)")
    bad = {k: v for k, v in ref2rows.items() if len(v) > 1}
    say(f"  reference loci spread over >1 row : {len(bad):,}")
    for k, v in list(bad.items())[:10]:
        say(f"      {k} -> {len(v)} rows")
    say("")
    say("[ 6 ] GENE UNIQUENESS")
    dup = sum(1 for v in gene_seen.values() if v > 1)
    say(f"  distinct genes                 : {len(gene_seen):,}")
    say(f"  genes appearing more than once : {dup:,}  (must be 0)")
    say("")
    if pos and multi_cells_genes:
        say("[ 7 ] SPAN PURITY   (multi-gene cells: are they tandem or dispersed?)")
        spans = []
        for gs in multi_cells_genes:
            co = [pos[g] for g in gs if g in pos]
            if len(co) < 2:
                continue
            s = [x[1] for x in co]
            spans.append(max(s) - min(s))
        if spans:
            sp = np.array(spans, dtype=float)
            for lo, hi, lab in ((0, 10_000, "< 10 kb   (tandem array)"),
                                (10_000, 100_000, "10-100 kb"),
                                (100_000, 1_000_000, "0.1-1 Mb"),
                                (1_000_000, 10**12, "> 1 Mb    (dispersed)")):
                m = ((sp >= lo) & (sp < hi)).sum()
                say(f"      {lab:26s}: {m:8,} ({100.0*m/len(sp):5.1f}%)")
            say(f"      median span = {np.median(sp):,.0f} bp")
        say("")
    say("=" * 88)
    say(" READING GUIDE")
    say("   * section 2 counts copies per CELL.  A high rate is good, but a")
    say("     multi-gene cell is NOT an error by itself: it means one subgenome")
    say("     carries several copies at that locus.  Section 3 shows whether the")
    say("     row is otherwise complete.")
    say("   * over-merge (section 5) must be empty for a clean reference frame.")
    say("   * genes appearing twice (section 6) must be 0: the matrix is a")
    say("     partition of the gene set.")
    say("   * synthetic rows have no reference orthologue.  They are positional")
    say("     annotations, not verified loci -- treat them separately downstream.")
    say("")
    say(" THREE THINGS CALLED 'SINGLE COPY' -- do not conflate them")
    say("   singleton row    one filled cell.  A COVERAGE statement: the other")
    say("                    subgenomes have no gene at this locus -- either lost,")
    say("                    or present but not placed.  Section 4 counts these.")
    say("   single-copy gene no paralogue anywhere in the genome (the usual meaning")
    say("                    of the term).  Implies a singleton row, not the reverse.")
    say("   complete 1:1     every column filled, one gene each: one copy per")
    say("                    subgenome and no extra duplication.  Sections 2 + 3.")
    say("=" * 88)

    if args.report:
        with open(args.report, "w") as fh:
            fh.write("\n".join(out) + "\n")
        print(f"\nreport written to {args.report}")


if __name__ == "__main__":
    main()
