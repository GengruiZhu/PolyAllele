#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The collision rule, and the surplus rule, pinned down.

Two genes on one chromosome can project to the same reference locus.  Sometimes
that is a real tandem array; sometimes it is a projection artefact, because
intervals in reference-locus space can overlap when the projection is
non-monotonic.  The two are told apart by GENE ORDER.

The first version of this test only PRINTED the outcome instead of asserting it,
which is how a behaviour change went unnoticed when `distribute_surplus` was later
flipped to false.  A test that cannot fail is not a test; every case below now has
an expected value.

  CASE A  anchors on ADJACENT reference loci: no room at all, so every gene
          between them gets its own slot
  CASE B  more genes than loci, distribute_surplus = false (the default): the
          first g_gap genes take the loci in order, the surplus gets slots
  CASE C  the same interval with distribute_surplus = true: the surplus is spread
          over the interval's loci and neighbours DO share a cell
  CASE D  fewer genes than loci: proportional spacing, one per locus, no collision
"""
import os
import sys
from collections import Counter

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

FAILED = []


def check(label, got, want):
    ok = got == want
    print(f"  {'ok  ' if ok else 'FAIL'} {label}")
    if not ok:
        print(f"         got  {got!r}")
        print(f"         want {want!r}")
        FAILED.append(label)


def run(n_loci, n_genes, anchors, window=5, distribute=False):
    eru_list = {"Chr01": ["g%d" % i for i in range(n_loci)]}
    eru_order = {g: ("Chr01", i) for i, g in enumerate(eru_list["Chr01"])}
    genes = ["A"] + ["x%d" % i for i in range(n_genes)] + ["B"]
    tier1 = {"A": (anchors[0], 1.0), "B": (anchors[1], 1.0)}
    res, _st, _unp = P.interpolate("ChrX", genes, tier1, eru_order, eru_list,
                                   50, 50, True, window, distribute)
    by_locus = Counter(res[g][0] for g in genes if g in res)
    multi = sum(v for v in by_locus.values() if v > 1)
    slots = sum(1 for g in res
                if "|iv|" in res[g][0] or "|same|" in res[g][0])
    return multi, slots


def main():
    print("=" * 78)
    print(" COLLISION AND SURPLUS RULES")
    print("=" * 78)

    print("\n CASE A -- anchors on adjacent reference loci (g_gap = 0)")
    check("no collision, every surplus gene gets a slot",
          run(7, 4, ("g2", "g3")), (0, 4))

    print("\n CASE B -- 8 genes over 5 loci, distribute_surplus = false")
    check("5 take the loci, 3 get slots, no collision",
          run(7, 8, ("g0", "g6")), (0, 3))

    print("\n CASE C -- the same interval, distribute_surplus = true")
    multi, slots = run(7, 8, ("g0", "g6"), distribute=True)
    check("the surplus is absorbed, so nothing is pushed to a slot", slots, 0)
    check("and neighbours end up sharing a cell", multi > 0, True)

    print("\n CASE D -- 5 genes over 18 loci")
    check("proportional spacing, one per locus, no slot",
          run(20, 5, ("g0", "g19")), (0, 0))

    print("\n" + "=" * 78)
    print("PASS" if not FAILED else f"FAIL ({len(FAILED)}): {FAILED}")
    print("=" * 78)
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
