#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Regression test for the v4.2d crash:

    KeyError: 'ROC_gene232936'   in stage_repair

A gene that entered NO row at all is not in `gene_meta` -- it never got one.
`stage_repair` nonetheless offered those genes up for repair (they are the
"interval too wide" ones, and repairing them is the whole point), then wrote
`gene_meta[h][1]` for a key that did not exist.

This drives `stage_repair` with a mocked aligner and asserts that a gene coming
from the unplaced set is moved onto its allelic locus, that it gains a gene_meta
entry, and that it leaves the unplaced set.
"""
import csv
import os
import shutil
import sys
import tempfile

def _find_bin():
    """Locate bin/ wherever this file has been copied to.

    Hard-coding the package directory name made every test break the moment the
    package was renamed -- which is exactly what a release does.  Look in every
    place the file can legitimately be: inside the package (test/), or beside it
    in the development tree.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (os.path.join(here, "..", "bin"),
                 os.path.join(here, "bin"),
                 os.path.join(here, "PolyAllele_v4", "bin"),
                 os.path.join(here, "..", "PolyAllele_v4", "bin")):
        if os.path.exists(os.path.join(cand, "polyallele.py")):
            return os.path.abspath(cand)
    raise SystemExit("cannot locate bin/polyallele.py relative to " + here)


sys.path.insert(0, _find_bin())
import polyallele as P  # noqa: E402


class Cfg:
    def __init__(self, wd):
        self.wd = wd
        self.od = wd
        self._s = {"_eru_order": {}}

    def s(self, sec, key, dv=""):
        return {"evalue": "1e-5"}.get(key, dv)

    def i(self, sec, key, dv=0):
        return {"threads": 2, "max_target_seqs": 20,
                "max_target_seqs_repair": 50, "shards": 1}.get(key, dv)

    def f(self, sec, key, dv=0.0):
        return {"min_identity": 90.0}.get(key, dv)

    def b(self, sec, key, dv=False):
        return dv

    def p(self, sec, key, dv=""):
        return dv

    def get(self, k, d=None):
        return self._s.get(k, d)

    def set(self, k, v):
        self._s[k] = v


def main():
    d = tempfile.mkdtemp()
    os.makedirs(os.path.join(d, "target"))
    os.makedirs(os.path.join(d, "reference"))
    os.makedirs(os.path.join(d, "tmp"))

    # two chromosomes of one group, plus a reference
    with open(os.path.join(d, "target", "Chr1A.bed"), "w") as o:
        o.write("Chr1A\t100\t200\tgA\n")
    with open(os.path.join(d, "target", "Chr1B.bed"), "w") as o:
        o.write("Chr1B\t100\t200\tgU\n")          # the unplaced gene
    for c in ("Chr1A", "Chr1B"):
        with open(os.path.join(d, "target", c + ".pep"), "w") as o:
            o.write(">%s\nMKVAAFLGLLSTTQMKVAAFLGLLSTTQMKVAAFLGLLSTTQMKVAAFLGLLSTTQMKVAAFLGLLSTTQMKVAAFLGLLSTTQMKVAAFLGLLSTTQMKVAA\n"
                    % ("gA" if c == "Chr1A" else "gU"))
    with open(os.path.join(d, "reference", "Chr01.bed"), "w") as o:
        o.write("Chr01\t100\t200\tL1\nChr01\t300\t400\tL2\n")
    with open(os.path.join(d, "manifest.tsv"), "w") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["Group", "RefChrom", "Chromosome", "CDS", "BED", "PEP"])
        for c in ("Chr1A", "Chr1B"):
            w.writerow(["G1", "Chr01", c, "", "", os.path.join(d, "target", c + ".pep")])

    cfg = Cfg(d)
    eru = {"L1": ("Chr01", 0), "L2": ("Chr01", 1)}
    cfg.set("_eru_order", eru)
    # gA sits at L1; gU is UNPLACED but bracketed by anchors at loci 0 and 1
    cfg.set("_gene_loc", {"gA": "L1"})
    cfg.set("_gene_meta", {"gA": ("anchor", "", "G1", "Chr1A", False, None)})
    cfg.set("_unplaced", {"gU": ("Chr01|wide|0-1", (0, 1),
                                 "interval_wider_than_max_interval", "G1", "Chr1B")})

    # mocked aligner: gU's only hit is gA at 97% over 28 aa
    P.which = lambda prog: "/usr/bin/" + prog
    P._run_besthit = lambda *a, **k: True
    P._read_hits = lambda path: {"gU": [("gA", 97.0, 108, 200.0)]}

    P.stage_repair(cfg)

    loc = cfg.get("_gene_loc")
    meta = cfg.get("_gene_meta")
    unpl = cfg.get("_unplaced")
    print()
    print("gU locus        :", loc.get("gU"))
    print("gU in gene_meta :", "gU" in meta)
    print("gU meta         :", meta.get("gU"))
    print("gU still unplaced:", "gU" in unpl)

    ok = (loc.get("gU") == "L1" and "gU" in meta
          and not meta["gU"][4]                       # not synthetic any more
          and "gU" not in unpl)
    print()
    print("PASS" if ok else "FAIL")

    # a second call must be idempotent (nothing left to repair)
    P.stage_repair(cfg)
    print("second pass was a no-op:",
          cfg.get("_gene_loc").get("gU") == "L1")

    # ------------------------------------------------------------------ gate 2
    # A move must not double up a cell on the gene's OWN chromosome.
    # gA (Chr1A) sits at L1; gC (Chr1B) also sits at L1, so cell (L1, Chr1B) is
    # taken; gV is on Chr1B and its only >=90% partner is gA on Chr1A -- a
    # DIFFERENT chromosome, so the paralog gate does not fire and only the
    # target-cell gate can stop it.
    print()
    print("--- target-cell gate ---")
    cfg2 = Cfg(d)
    cfg2.set("_eru_order", eru)
    cfg2.set("_gene_loc", {"gA": "L1", "gC": "L1"})
    cfg2.set("_gene_meta", {"gA": ("anchor", "", "G1", "Chr1A", False, None),
                            "gC": ("interp_exact", "", "G1", "Chr1B", False, (0, 1))})
    cfg2.set("_unplaced", {"gV": ("Chr01|wide|0-1", (0, 1),
                                  "interval_wider_than_max_interval", "G1", "Chr1B")})
    P._read_hits = lambda path: {"gV": [("gA", 97.0, 108, 200.0)]}
    P.stage_repair(cfg2)
    moved = cfg2.get("_gene_loc").get("gV")
    print()
    print("gV locus            :", moved)
    print("gV still unplaced   :", "gV" in cfg2.get("_unplaced"))
    ok2 = moved is None and "gV" in cfg2.get("_unplaced")
    print("PASS" if ok2 else "FAIL")

    shutil.rmtree(d)


if __name__ == "__main__":
    main()
