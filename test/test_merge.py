#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stage `merge` -- every gate, one case each.

A row holding one gene is read downstream as "nothing else in the genome is this
gene's allele".  That is only true when the gene really is alone.  When its
homoeolog was placed in ANOTHER row, the row lies, and the gene has to move.

The move is refused unless ALL of these hold, and each refusal has its own case:

  CASE 1  merge          alone, partner in the table, different chromosome,
                         destination cell free
  CASE 2  refuse         destination cell already holds a gene of the mover's
                         own chromosome -- filling it would put two genes of one
                         chromosome in one cell, the over-merge this pipeline
                         exists to avoid
  CASE 3  refuse         partner sits on the mover's OWN chromosome: a paralog,
                         not an allele
  CASE 4  refuse         different homology group
  CASE 5  refuse         partner is not in the table at all
  CASE 6  untouched      a gene that was never alone is never moved
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


class Cfg:
    def __init__(self, wd):
        self.wd = wd
        self.od = wd
        self._s = {}

    def s(self, sec, key, dv=""):
        return {"evalue": "1e-5"}.get(key, dv)

    def i(self, sec, key, dv=0):
        return {"threads": 2, "max_target_seqs": 20,
                "max_target_seqs_merge": 50, "shards": 1}.get(key, dv)

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
    os.makedirs(os.path.join(d, "tmp"), exist_ok=True)

    cfg = Cfg(d)
    #  gene : (source, interval, group, chromosome, synthetic, bounds)
    meta = {
        "gA": ("anchor", "", "G1", "Chr1A", False, None),
        "gZ": ("anchor", "", "G1", "Chr1C", False, None),
        "gQ": ("anchor", "", "G1", "Chr1B", False, None),
        "gX": ("interp_synthetic", "iv", "G1", "Chr1B", True, None),
        "gY": ("interp_synthetic", "iv", "G1", "Chr1C", True, None),
        "gW": ("interp_synthetic", "iv", "G1", "Chr1B", True, None),
        "gV": ("interp_synthetic", "iv", "G1", "Chr1D", True, None),
        "gU": ("interp_synthetic", "iv", "G2", "Chr2A", True, None),
        "gM": ("anchor", "", "G1", "Chr1E", False, None),
        "gN": ("anchor", "", "G1", "Chr1F", False, None),
        # CASE 7 -- a tier-1 anchor, alone in its row, with a perfect partner
        # available and a free cell.  Position outranks similarity here.
        "gR": ("anchor", "", "G1", "Chr1I", False, None),
    }
    #  L1 is a two-cell row, so it can receive; every other row holds one gene
    loc = {"gA": "L1", "gZ": "L1",
           "gQ": "L6",
           "gX": "L2", "gY": "L3", "gW": "L4", "gV": "L5", "gU": "L7",
           "gM": "L8", "gN": "L8", "gR": "L9"}
    cfg.set("_gene_loc", dict(loc))
    cfg.set("_gene_meta", dict(meta))
    cfg.set("_unplaced", {})

    hits = {
        # CASE 1 -- partner gA on Chr1A, mover on Chr1B, cell (L1, Chr1B) free
        "gX": [("gA", 96.0, 120, 300.0)],
        # CASE 2 -- partner gA, but cell (L1, Chr1C) already holds gZ
        "gY": [("gA", 96.0, 120, 300.0)],
        # CASE 3 -- partner gQ sits on the mover's own chromosome (Chr1B)
        "gW": [("gQ", 99.0, 120, 350.0)],
        # CASE 4 -- partner gA is in another homology group
        "gU": [("gA", 97.0, 120, 320.0)],
        # CASE 5 -- the best hit is not in the table
        "gV": [("gNOPE", 99.0, 120, 400.0)],
        # CASE 7 -- a perfect, fully eligible partner, and still no move
        "gR": [("gA", 100.0, 120, 500.0)],
    }

    P.which = lambda prog: "/usr/bin/" + prog
    P.load_manifest = lambda wd: [("G1", "refA", c, c + ".cds", c + ".pep")
                                  for c in ("Chr1A", "Chr1B", "Chr1C", "Chr1D",
                                            "Chr1E", "Chr1F", "Chr2A")]
    P._write_pep = lambda pep_of, genes, path: len(list(genes))
    P._run_besthit = lambda *a, **k: True
    P._read_hits = lambda path: hits

    print("=" * 78)
    print(" STAGE merge")
    print("=" * 78)
    P.stage_merge(cfg)

    gl = cfg.get("_gene_loc")
    gm = cfg.get("_gene_meta")
    merged = cfg.get("_merged")
    print()
    print("--- assertions ---")
    check("CASE 1 the lone gene moved to the partner's row", gl.get("gX"), "L1")
    check("CASE 1 it is recorded as merged",
          gm["gX"][0], "merged_homoeolog")
    check("CASE 1 its chromosome did not change", gm["gX"][3], "Chr1B")
    check("CASE 1 its own row is gone",
          [g for g, L in gl.items() if L == "L2"], [])
    check("CASE 2 cell taken -> stayed put", gl.get("gY"), "L3")
    check("CASE 3 same-chromosome partner -> stayed put", gl.get("gW"), "L4")
    check("CASE 4 different group -> stayed put", gl.get("gU"), "L7")
    check("CASE 5 partner not in the table -> stayed put", gl.get("gV"), "L5")
    check("CASE 6 a gene that was never alone never moved", gl.get("gA"), "L1")
    check("CASE 6 ... and neither did its neighbour", gl.get("gM"), "L8")
    check("CASE 7 an anchored gene stays where collinearity put it",
          gl.get("gR"), "L9")
    check("CASE 7 ... and keeps its anchor source", gm["gR"][0], "anchor")
    check("merged set is exactly the one gene", sorted(merged), ["gX"])
    check("L1 now holds three genes",
          sorted(g for g, L in gl.items() if L == "L1"), ["gA", "gX", "gZ"])

    print()
    print("--- idempotence: a second call must move nothing ---")
    before = dict(gl)
    P.stage_merge(cfg)
    check("second call changed nothing", gl == before, True)

    shutil.rmtree(d, ignore_errors=True)
    print()
    print("=" * 78)
    print("PASS" if not FAILED else f"FAIL ({len(FAILED)}): {FAILED}")
    print("=" * 78)
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
