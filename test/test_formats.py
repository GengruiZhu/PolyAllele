#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Input-format compatibility.

A gene ID must match between the target BED, the target sequences and the JCVI
anchors, or NOTHING is anchored -- and the failure surfaces much later as an empty
result, not as an error.  Annotation pipelines disagree about the exact spelling,
so these tests pin down what the readers accept:

  * gene IDs written `gene`, `gene-P`, `gene-T`, `"gene-T"`  -> one gene
  * GFF3 with `gene` features and `ID=`
  * GTF with only `transcript` features and `gene_id "..."` (quoted)
  * a GFF holding only `CDS` records
  * CRLF line endings
  * BED separated by spaces, with a `track` line, with extra columns
  * a legacy FASTA comment line
"""
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "bin"))
import polyallele as P  # noqa: E402

FAILED = []


def check(label, got, want):
    ok = got == want
    print(f"  {'ok  ' if ok else 'FAIL'} {label}")
    if not ok:
        print(f"         got  {got!r}")
        print(f"         want {want!r}")
        FAILED.append(label)


def main():
    d = tempfile.mkdtemp()
    try:
        print("=" * 78)
        print(" ID NORMALISATION")
        print("=" * 78)
        P.set_id_strip(r"-[A-Za-z]$")
        check("plain id unchanged", P.norm_id("ROC_gene123"), "ROC_gene123")
        check("-P stripped", P.norm_id("ROC_gene123-P"), "ROC_gene123")
        check("-T stripped", P.norm_id("ROC_gene123-T"), "ROC_gene123")
        check("quotes stripped too", P.norm_id('"ROC_gene9-T"'), "ROC_gene9")
        check("missing value tolerated", P.norm_id("  g1  "), "g1")
        P.set_id_strip("")
        check("pattern empty -> exact match", P.norm_id("g1-P"), "g1-P")
        P.set_id_strip(r"-[A-Za-z]$")

        print()
        print("=" * 78)
        print(" GFF / GTF")
        print("=" * 78)
        g = os.path.join(d, "a.gff")
        with open(g, "w", newline="") as o:
            o.write("##gff-version 3\r\n")
            o.write("chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=geneA;Name=x\r\n")
            o.write("chr1\tsrc\tmRNA\t100\t200\t.\t+\t.\tID=geneA.1;Parent=geneA\r\n")
        check("GFF3 gene + ID, CRLF",
              P.parse_gff_genes(g), {"geneA": ("chr1", 100, 200)})

        g = os.path.join(d, "b.gtf")
        with open(g, "w") as o:
            o.write('chr1\tsrc\ttranscript\t5\t90\t.\t-\t.\t'
                    'gene_id "gB"; transcript_id "gB.1";\n')
        check("GTF transcript + quoted gene_id",
              P.parse_gff_genes(g), {"gB": ("chr1", 5, 90)})

        g = os.path.join(d, "c.gff")
        with open(g, "w") as o:
            o.write("chr1\tsrc\tCDS\t7\t60\t.\t+\t0\tID=cdsC\n")
        check("GFF holding only CDS",
              P.parse_gff_genes(g), {"cdsC": ("chr1", 7, 60)})

        g = os.path.join(d, "d.gff")
        with open(g, "w") as o:
            o.write("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tName=geneD\n")
        check("Name= instead of ID=",
              P.parse_gff_genes(g), {"geneD": ("chr1", 1, 9)})

        g = os.path.join(d, "e.gff")
        with open(g, "w") as o:
            o.write("# only comments\n")
        check("no genes -> empty, not a crash", P.parse_gff_genes(g), {})

        print()
        print("=" * 78)
        print(" BED")
        print("=" * 78)
        b = os.path.join(d, "x.bed")
        with open(b, "w") as o:
            o.write("track name=x\n")
            o.write("chr1 10 20 gD-P extra\n")
            o.write("chr1 30 40 gE\n")
        check("space separated, track line, extra columns, -P stripped",
              P.read_bed(b), [(10, "gD"), (30, "gE")])

        b = os.path.join(d, "y.bed")
        with open(b, "w", newline="") as o:
            o.write("chr1\t30\t40\tgE\r\n")
            o.write("chr1\t10\t20\tgD\r\n")
        check("tab separated with CRLF, sorted by position",
              P.read_bed(b), [(10, "gD"), (30, "gE")])

        print()
        print("=" * 78)
        print(" FASTA / ANCHORS")
        print("=" * 78)
        f = os.path.join(d, "z.pep")
        with open(f, "w", newline="") as o:
            o.write("; a legacy comment\r\n")
            o.write(">gD-P some description here\r\nMKV\r\n\r\nAA\r\n")
            o.write(">gE\nTTT\n")
        check("legacy comment, CRLF, blank line, description trimmed",
              list(P.read_fasta(f)), [("gD", "MKVAA"), ("gE", "TTT")])

        a = os.path.join(d, "x.anchors")
        with open(a, "w") as o:
            o.write("# comment\n")
            o.write("gD-P\trefA\t0.9\n")
            o.write("gE\trefB\n")
        check("anchors: comment, -P stripped, missing score tolerated",
              P.read_anchors(a), [("gD", "refA", 0.9), ("gE", "refB", 0.0)])

        # ------------------------------------------------------------------ LF
        # csv.writer defaults to CRLF.  A CRLF TSV breaks `awk '{print $NF}'`,
        # `cut`, and any `line.split("\t")[-1]` on Linux -- the last field keeps a
        # stray \r.  It broke this project's own profiler, which read a column
        # header as 'Source\r' and reported every value as missing.  Every writer
        # must therefore say lineterminator="\n", and this test says so.
        print()
        print("=" * 78)
        print(" OUTPUT LINE ENDINGS")
        print("=" * 78)
        bindir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin")
        offenders = []
        for fn in sorted(os.listdir(bindir)):
            if not fn.endswith(".py"):
                continue
            text = open(os.path.join(bindir, fn), encoding="utf-8").read()
            for i, line in enumerate(text.splitlines(), 1):
                if "csv.writer(" in line and "lineterminator" not in line:
                    offenders.append(f"{fn}:{i}")
        check("every csv.writer sets lineterminator", offenders, [])

        import csv as _csv
        import io as _io
        buf = _io.StringIO()
        _csv.writer(buf, delimiter="\t", lineterminator="\n").writerow(["a", "b"])
        check("and the result really is LF, not CRLF",
              buf.getvalue(), "a\tb\n")
    finally:
        shutil.rmtree(d, ignore_errors=True)

    print()
    print("=" * 78)
    print("PASS" if not FAILED else f"FAIL ({len(FAILED)}): {FAILED}")
    print("=" * 78)
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
