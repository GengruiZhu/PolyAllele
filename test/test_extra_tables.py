#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test the tables that sit beside the matrix.

Five behaviours, one case each:

  CASE 1  an unplaced gene whose relative is in the table joins that locus
  CASE 2  two unplaced genes that are each other's relative share a NEW locus,
          so neither of them looks single copy
  CASE 3  an unplaced gene with no relative gets its own locus and IS single copy
  CASE 4  RowIsSingleCopy needs BOTH conditions: the locus holds one gene AND
          that gene has no relative.  A gene alone at its locus but with a
          relative elsewhere must NOT be flagged -- this is the distinction the
          whole file exists to make
  CASE 5  competition clusters use their own threshold and ignore position, and
          every output file is written with LF, not CRLF
"""
import csv
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "bin", "extra_tables.py")
FAILED = []


def check(label, got, want):
    ok = got == want
    print(f"  {'ok  ' if ok else 'FAIL'} {label}")
    if not ok:
        print(f"         got  {got!r}")
        print(f"         want {want!r}")
        FAILED.append(label)


def write(path, text):
    with open(path, "w", newline="") as o:
        o.write(text)


def read_tsv(path):
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t")
                if r and not r[0].startswith("#")]
    head, body = rows[0], rows[1:]
    return head, [dict(zip(head, r)) for r in body]


def main():
    d = tempfile.mkdtemp()
    try:
        m = os.path.join(d, "m.tsv")
        write(m, "LocusID\tHomologyGroup\tAllele_A\tAllele_B\n"
                 "L1\tG1\tgA1\tgB1\n"
                 "L2\tG1\tNA\tgB2\n")
        lng = os.path.join(d, "l.tsv")
        write(lng, "LocusID\tHomologyGroup\tChromosome\tGeneID\tSource\t"
                   "Subgenome\tRank\n"
                   "L1\tG1\tChr1A\tgA1\tanchor\tA\t1\n"
                   "L1\tG1\tChr1B\tgB1\tanchor\tB\t1\n"
                   "L2\tG1\tChr1B\tgB2\tanchor\tB\t1\n")
        unp = os.path.join(d, "u.tsv")
        write(unp, "GeneID\tChromosome\tGeneOrderInChrom\tStart\tEnd\n"
                   "gU1\tChr1C\t0\t1\t2\n"
                   "gU2\tChr1D\t0\t1\t2\n"
                   "gU3\tChr1E\t0\t1\t2\n"
                   "gU4\tChr1F\t0\t1\t2\n"
                   "gU5\tChr1G\t0\t1\t2\n")
        av = os.path.join(d, "a.tsv")
        write(av, "gU1\tgA1\t98.0\t100\n"
                  "gU2\tgU3\t97.0\t100\n"
                  "gU3\tgU2\t97.0\t100\n"
                  "gA1\tgB1\t99.0\t100\n"
                  "gA1\tgU1\t98.0\t100\n"
                  "gB2\tgB2\t100.0\t100\n"
                  # gU4 has NO hit at all -> single copy.
                  # gU5 has only a weak one, and that still counts: CopyClass
                  # asks whether ANYTHING in the proteome resembles the gene, at
                  # any identity BLASTP reports.  The threshold for "related" is
                  # the aligner's, not this script's -- which is exactly why the
                  # table also carries BestAnyPident, so a stricter reading can
                  # be applied afterwards.
                  "gU5\tgZ1\t40.0\t100\n")

        r = subprocess.run([sys.executable, SCRIPT,
                            "--matrix", m, "--long", lng, "--unplaced", unp,
                            "--allvall", av, "--out-dir", d, "--prefix", "t",
                            "--identity", "90", "--competition-identity", "95",
                            "--min-alignment", "60"],
                           capture_output=True, text=True)
        print("=" * 78)
        print(" extra_tables.py")
        print("=" * 78)
        if r.returncode != 0:
            print(r.stdout)
            print(r.stderr)
            check("script ran", r.returncode, 0)
            return 1

        head, g2c = read_tsv(os.path.join(d, "t.gene2copyclass.tsv"))
        by = {x["GeneID"]: x for x in g2c}
        check("CASE 1 the unplaced gene joined the locus of its relative",
              by["gU1"]["LocusID"], "L1")
        check("CASE 1 the locus now holds three genes",
              by["gA1"]["LocusSize"], "3")
        check("CASE 2 two homeless genes share a new locus",
              by["gU2"]["LocusID"] == by["gU3"]["LocusID"], True)
        check("CASE 2 ... and it holds exactly those two",
              by["gU2"]["LocusSize"], "2")
        check("CASE 2 neither of them is single copy",
              (by["gU2"]["CopyClass"], by["gU3"]["CopyClass"]),
              ("multi_copy", "multi_copy"))
        check("CASE 3 a gene with no relative at all",
              by["gU4"]["CopyClass"], "single_copy")
        check("CASE 3 ... gets its own locus", by["gU4"]["LocusSize"], "1")
        check("CASE 3 ... and is flagged as a clean single-copy row",
              by["gU4"]["RowIsSingleCopy"], "1")
        check("CASE 4 alone at its locus, no relative -> flagged",
              by["gB2"]["RowIsSingleCopy"], "1")
        check("CASE 4 alone at its locus BUT has a relative -> not flagged",
              (by["gU1"]["LocusSize"], by["gU1"]["CopyClass"],
               by["gU1"]["RowIsSingleCopy"]),
              ("3", "multi_copy", "0"))
        check("CASE 4 ... the flag and the class are independent",
              by["gB1"]["RowIsSingleCopy"], "0")
        check("CASE 5 competition ignores position: gU1 joins L1's genes",
              by["gU1"]["CompetitionClusterID"], by["gA1"]["CompetitionClusterID"])
        check("CASE 5 ... and the pair forms its own competition cluster",
              by["gU2"]["CompetitionClusterID"] ==
              by["gU3"]["CompetitionClusterID"]
              and by["gU2"]["CompetitionClusterID"] !=
              by["gA1"]["CompetitionClusterID"], True)
        check("CASE 5 a gene nothing can be confused with stands alone",
              by["gB2"]["CompetitionClusterSize"], "1")

        head, l2g = read_tsv(os.path.join(d, "t.locus2gene.tsv"))
        check("locus2gene covers every gene",
              sorted(x["GeneID"] for x in l2g),
              ["gA1", "gB1", "gB2", "gU1", "gU2", "gU3", "gU4", "gU5"])
        head, c2g = read_tsv(os.path.join(d, "t.competition2gene.tsv"))
        check("competition2gene covers every gene too",
              len(c2g), 8)

        for f in ("t.gene2copyclass.tsv", "t.locus2gene.tsv",
                  "t.competition2gene.tsv", "t.unplaced_matrix.tsv"):
            with open(os.path.join(d, f), "rb") as fh:
                blob = fh.read()
            check(f"{f} uses LF", b"\r\n" in blob, False)
    finally:
        shutil.rmtree(d, ignore_errors=True)

    print()
    print("=" * 78)
    print("PASS" if not FAILED else f"FAIL ({len(FAILED)}): {FAILED}")
    print("=" * 78)
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
