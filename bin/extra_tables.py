#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PolyAllele v4.4 -- the tables that sit NEXT TO the allele matrix.

The matrix answers "which ancestral locus does this gene belong to".  It cannot
answer "is this gene single copy", and the two questions have to be kept apart:

  * a gene can be alone AT ITS LOCUS and still have a >=90% relative AT ANOTHER
    locus.  Such a gene is not single copy, and its row is nonetheless correct;
  * the downstream quantification tool does not read these labels.  It reads the
    rows of a cluster file and treats a one-gene row as single copy, so two genes
    that are really brothers but sit in separate one-gene rows will NOT compete,
    and their expression comes out inflated.

Those two facts pull in opposite directions, so this script emits BOTH groupings
rather than pretending one can serve both:

  LOCUS cluster        an ancestral locus.  Position-limited: an unplaced gene is
                       attached to the row holding its best relative, and a cell
                       that is already occupied is not a reason to refuse --
                       this table counts copies, it does not claim allelism.  A
                       locus holds tens of genes, never tens of thousands.

  COMPETITION cluster  genes whose reads cannot be told apart: connected
                       components of the >=T identity graph, position ignored.
                       T defaults to 95, because at 150 bp, 90% identity leaves
                       about fifteen mismatches and an aligner resolves those;
                       the ambiguity that matters is higher up.  THIS is the
                       grouping to hand to a quantification tool.

Files written, all with LF endings (csv.writer defaults to CRLF, which breaks
awk, cut, and any `line.split("\\t")[-1]` on Linux -- that is not theoretical, it
silently corrupted this project's own profiler):

  <prefix>.unplaced_matrix.tsv    genes in no row, arranged as a matrix
  <prefix>.locus2gene.tsv         LocusID <TAB> GeneID            (position)
  <prefix>.competition2gene.tsv   ClusterID <TAB> GeneID          (quantification)
  <prefix>.gene2copyclass.tsv     one line per gene, both groupings and labels
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict

NULLS = {"NA", "NAN", "", ".", "-"}
LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def log(m):
    print(m, flush=True)


def split_cell(v):
    v = str(v).strip()
    return [] if v.upper() in NULLS else [x.strip() for x in v.split(",") if x.strip()]


def base_id(g):
    return re.sub(r"-[A-Za-z]$", "", str(g).strip())


def group_of(chrom):
    m = re.match(r"^Chr(\d+)", chrom or "")
    return m.group(1) if m else (chrom or "?")


def rows_of(path):
    """(header, iterator of field lists), skipping `#` comment lines."""
    fh = open(path, newline="", errors="replace")
    head = None
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        head = line.rstrip("\n").split("\t")
        break

    def gen():
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            yield line.rstrip("\n").split("\t")
    return head, gen()


class Union:
    __slots__ = ("p",)

    def __init__(self):
        self.p = {}

    def find(self, x):
        p = self.p
        if x not in p:
            p[x] = x
            return x
        r = x
        while p[r] != r:
            r = p[r]
        while p[x] != r:
            p[x], x = r, p[x]
        return r

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description="PolyAllele extra tables")
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--long", required=True)
    ap.add_argument("--unplaced", required=True)
    ap.add_argument("--allvall", required=True,
                    help="genome-wide all-vs-all in BLAST outfmt 6, first four "
                         "columns qseqid sseqid pident length.  MUST come from "
                         "BLASTP: measured against BLASTP, DIAMOND --sensitive "
                         "called 24.58%% of its single_copy genes wrongly.")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default="polyallele")
    ap.add_argument("--identity", type=float, default=90.0,
                    help="identity for attaching an unplaced gene to a locus")
    ap.add_argument("--competition-identity", type=float, default=95.0,
                    help="identity above which two genes fight over one read")
    ap.add_argument("--min-alignment", type=int, default=60)
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    P = a.prefix

    log("=" * 92)
    log(" PolyAllele -- the tables beside the matrix")
    log("=" * 92)

    # ---- the table ---------------------------------------------------------
    head, gen = rows_of(a.long)
    ci = {c: i for i, c in enumerate(head)}
    g2loc, g2chrom, g2src = {}, {}, {}
    for f in gen:
        g = f[ci["GeneID"]]
        g2loc[g] = f[ci["LocusID"]]
        g2chrom[g] = f[ci["Chromosome"]]
        g2src[g] = f[ci["Source"]] if "Source" in ci else "?"
    log(f"  table genes     : {len(g2loc):,}   loci: {len(set(g2loc.values())):,}")

    head, gen = rows_of(a.unplaced)
    ui = {c: i for i, c in enumerate(head)}
    unplaced, unp_chrom = [], {}
    for f in gen:
        g = f[ui["GeneID"]]
        unplaced.append(g)
        unp_chrom[g] = f[ui["Chromosome"]] if "Chromosome" in ui else ""
    log(f"  unplaced genes  : {len(unplaced):,}")

    allg = list(g2loc) + unplaced

    # ---- chromosome -> matrix column ---------------------------------------
    mhead, mgen = rows_of(a.matrix)
    mi = {c: i for i, c in enumerate(mhead)}
    acols = [c for c in mhead if c.startswith("Allele_")]
    chrom_letter = {}
    for f in mgen:
        for c in acols:
            for g in split_cell(f[mi[c]]):
                ch = g2chrom.get(g)
                if ch:
                    chrom_letter.setdefault(ch, c)
    ncols = max((LETTERS.index(c[-1]) for c in chrom_letter.values()
                 if c[-1] in LETTERS), default=0) + 1
    cols = [f"Allele_{LETTERS[k]}" for k in range(ncols)]
    log(f"  chromosome columns : {ncols}")

    # ---- one streaming pass over the all-vs-all ----------------------------
    log(f"  reading {os.path.basename(a.allvall)} ...")
    G90 = Union()
    GC = Union()
    for g in allg:
        G90.find(g)
        GC.find(g)
    best_any = defaultdict(float)
    best_ge = defaultdict(float)
    partner = {}
    n90 = ncomp = 0
    with open(a.allvall, errors="replace") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or p[0] == p[1]:
                continue
            q, s = base_id(p[0]), base_id(p[1])
            try:
                pid, ln = float(p[2]), int(p[3])
            except ValueError:
                continue
            if ln < a.min_alignment:
                continue
            if pid > best_any[q]:
                best_any[q] = pid
            if pid >= a.competition_identity:
                GC.union(q, s)
                ncomp += 1
            if pid >= a.identity:
                G90.union(q, s)
                n90 += 1
                if pid > best_ge[q]:
                    best_ge[q] = pid
                    partner[q] = s
    log(f"    >= {a.identity:g}% edges : {n90:,}")
    log(f"    >= {a.competition_identity:g}% edges : {ncomp:,}")

    # ---- attach the unplaced genes to loci ---------------------------------
    row_of = dict(g2loc)
    for g in unplaced:
        row_of[g] = None
    attached = 0
    for g in unplaced:
        h = partner.get(g)
        if h and h in g2loc:
            row_of[g] = g2loc[h]
            attached += 1
    leftover = [g for g in unplaced if row_of[g] is None]
    uf = {g: g for g in leftover}

    def ufind(x):
        while uf[x] != x:
            uf[x] = uf[uf[x]]
            x = uf[x]
        return x

    for g in leftover:
        h = partner.get(g)
        if h and h in uf:
            ra, rb = ufind(g), ufind(h)
            if ra != rb:
                uf[rb] = ra
    for g in leftover:
        row_of[g] = f"NOANCHOR|{ufind(g)}"
    tmp = Counter(row_of[g] for g in leftover)
    paired = sum(1 for g in leftover if tmp[row_of[g]] > 1)
    alone_new = len(leftover) - paired
    locus_size = Counter(row_of[g] for g in allg)
    log(f"  attached to an existing locus : {attached:,}")
    log(f"  paired with another homeless  : {paired:,}")
    log(f"  alone in a new locus          : {alone_new:,}")

    # ---- competition groups -------------------------------------------------
    comp_of = {g: GC.find(g) for g in allg}
    comp_size = Counter(comp_of.values())

    # ---- 1. unplaced matrix -------------------------------------------------
    f_un = os.path.join(a.out_dir, f"{P}.unplaced_matrix.tsv")
    ucomp = defaultdict(list)
    for g in unplaced:
        ucomp[row_of[g]].append(g)
    tally = Counter()
    with open(f_un, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["# Rows here are LOCI.  An unplaced gene with no partner gets "
                    "its own row, and that row really is a single-copy gene."])
        w.writerow(["UnplacedLocusID", "HomologyGroup"] + cols
                   + ["NCopies", "Evidence", "Chromosomes"])
        for L in sorted(ucomp, key=lambda x: (-len(ucomp[x]), str(x))):
            gs = sorted(ucomp[L], key=lambda g: (unp_chrom.get(g, ""), g))
            row = {c: "NA" for c in cols}
            per = defaultdict(list)
            for g in gs:
                per[unp_chrom.get(g, "")].append(g)
            for chrom in sorted(per):
                col = chrom_letter.get(chrom)
                if col and row.get(col, "NA") == "NA":
                    row[col] = per[chrom][0]
            ev = ("similarity_pair"
                  if len(set(unp_chrom.get(g, "") for g in gs)) > 1
                  else "no_partner")
            tally[ev] += 1
            w.writerow([L, group_of(unp_chrom.get(gs[0], ""))]
                       + [row[c] for c in cols]
                       + [len(gs), ev,
                          ",".join(sorted(set(unp_chrom.get(g, "") for g in gs)))])
    log(f"  -> {os.path.basename(f_un)}   {sum(tally.values()):,} rows  "
        f"{dict(tally)}")

    # ---- 2. locus -> gene ---------------------------------------------------
    f_l2g = os.path.join(a.out_dir, f"{P}.locus2gene.tsv")
    with open(f_l2g, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["LocusID", "GeneID", "Chromosome", "Column", "Role",
                    "Source", "InTable"])
        seen = Counter()
        for g in allg:
            L = row_of[g]
            seen[L] += 1
            ch = g2chrom.get(g) or unp_chrom.get(g, "")
            w.writerow([L, g, ch, chrom_letter.get(ch, "?"),
                        "representative" if seen[L] == 1 else "extra_copy",
                        g2src.get(g, "unplaced"), int(g in g2loc)])
    log(f"  -> {os.path.basename(f_l2g)}   {len(allg):,} rows")

    # ---- 3. competition -> gene ---------------------------------------------
    f_c2g = os.path.join(a.out_dir, f"{P}.competition2gene.tsv")
    with open(f_c2g, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow([f"# Competition clusters: connected components of the "
                    f">= {a.competition_identity:g}% identity graph, position "
                    f"ignored."])
        w.writerow(["# Feed THIS file to a quantification tool.  A one-gene row "
                    "is a gene nothing can be confused with; brothers share a "
                    "row, so they compete instead of both counting in full."])
        w.writerow(["ClusterID", "GeneID", "ClusterSize", "InTable"])
        for g in sorted(allg, key=lambda x: (comp_of[x], x)):
            w.writerow([f"CC|{comp_of[g]}", g, comp_size[comp_of[g]],
                        int(g in g2loc)])
    log(f"  -> {os.path.basename(f_c2g)}   {len(allg):,} rows, "
        f"{len(comp_size):,} clusters, largest {max(comp_size.values()):,}")

    # ---- 4. per gene --------------------------------------------------------
    f_g2c = os.path.join(a.out_dir, f"{P}.gene2copyclass.tsv")
    single = 0
    clean_rows = 0
    with open(f_g2c, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["# CopyClass: does ANYTHING in the proteome resemble this "
                    "gene at all (BLASTP, any identity)?  single_copy means no."])
        w.writerow(["# LocusID / LocusSize : the ancestral locus, a POSITION.  "
                    "RowIsSingleCopy is the flag to read: 1 only when the locus "
                    "holds one gene AND that gene has no relative at all."])
        w.writerow([f"# CompetitionClusterID : genes whose reads cannot be told "
                    f"apart (>= {a.competition_identity:g}%).  Use THIS, not the "
                    f"locus, for quantification."])
        w.writerow(["GeneID", "CopyClass", "LocusID", "LocusSize",
                    "RowIsSingleCopy", "BestAnyPident", "BestNonSelfPident",
                    "CompetitionClusterID", "CompetitionClusterSize", "InTable"])
        for g in sorted(allg):
            b = best_any.get(g, 0.0)
            cls = "single_copy" if b == 0.0 else "multi_copy"
            single += cls == "single_copy"
            L = row_of[g]
            lsz = locus_size[L]
            cid = comp_of[g]
            flag = int(lsz == 1 and b == 0.0)
            clean_rows += flag
            w.writerow([g, cls, L, lsz, flag, f"{b:.1f}",
                        f"{best_ge.get(g, 0.0):.1f}", f"CC|{cid}",
                        comp_size[cid], int(g in g2loc)])

    log("")
    log("=" * 92)
    log(f"  genes                 : {len(allg):,}")
    log(f"  single_copy           : {single:,}  "
        f"({100.0*single/max(len(allg),1):.2f}%)   no hit at all")
    log(f"  loci                  : {len(locus_size):,}   "
        f"largest {max(locus_size.values()):,} genes")
    log(f"  loci holding 1 gene   : {sum(1 for v in locus_size.values() if v == 1):,}"
        f"   of which truly single copy: {clean_rows:,}")
    log(f"  competition clusters  : {len(comp_size):,}   "
        f"largest {max(comp_size.values()):,} genes")
    log("=" * 92)


if __name__ == "__main__":
    sys.exit(main())
