#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Local test of the query-sharding orchestration, with a fake aligner so it does
not need blastp.  Verifies: every query appears exactly once in the merged
output, the shard files are cleaned up, and the cache short-circuits a re-run."""
import os
import shutil
import sys
import tempfile
import time

def _find_bin():
    """Locate bin/ wherever this file has been copied to.

    Hard-coding the package directory name made every test break the moment the
    package was renamed -- which is exactly what a release does.  Look in the
    places the file can legitimately be: inside the package (test/), or beside it
    in the development tree.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (os.path.join(here, "..", "bin"),
                 os.path.join(here, "bin"),
                 os.path.join(here, "PolyAllele_v4", "bin"),
                 os.path.join(here, "..", "PolyAllele_v4", "bin"),
                 os.path.join(here, "..", "..", "PolyAllele_v4", "bin")):
        if os.path.exists(os.path.join(cand, "polyallele.py")):
            return os.path.abspath(cand)
    raise SystemExit("cannot locate bin/polyallele.py relative to " + here)


sys.path.insert(0, _find_bin())
import polyallele as P  # noqa: E402

FAKE = (
    "import sys\n"
    "o = open(sys.argv[2], 'w')\n"
    "for line in open(sys.argv[1]):\n"
    "    if line.startswith('>'):\n"
    "        o.write(line[1:].strip() + '\\t' + sys.argv[3] + '\\n')\n"
    "o.close()\n"
)


class Cfg:
    def s(self, sec, key, dv=""):
        return {"evalue": "1e-5"}.get(key, dv)

    def i(self, sec, key, dv=0):
        return {"threads": 12, "max_target_seqs": 50, "shards": 3}.get(key, dv)

    def b(self, sec, key, dv=False):
        return dv


def main():
    d = tempfile.mkdtemp()
    q = os.path.join(d, "q.pep")
    # Records must look like real proteins: a short header and a long sequence.
    # With tiny sequences the headers dominate the byte count, and a splitter that
    # only counts header bytes passes the test while failing on real data.
    with open(q, "w") as o:
        for i in range(300):
            o.write(">g%03d\n" % i)
            for _ in range(12):
                o.write("MKVAAFLGLLSTTQ\n")
    t = os.path.join(d, "t.pep")
    with open(t, "w") as o:
        o.write(">x\nMKV\n")

    # the split itself must be even
    shards = P._split_fasta(q, 4, d, "splittest")
    sizes = [os.path.getsize(s) for s in shards]
    total = os.path.getsize(q)
    even = all(s > 0 for s in sizes) and max(sizes) < 1.6 * (total / 4)
    print("split sizes       :", sizes, "of", total)
    print("split is even     :", even)
    for s in shards:
        os.remove(s)

    def fake_search(cfg, used, qf, db, out, threads, label=""):
        # one output line per query, tagged with the thread count it was given
        return [sys.executable, "-c", FAKE, qf, out, str(threads)]

    P._search_cmd = fake_search
    P._besthit_cmds = lambda cfg, tf, db, qf, o: (
        [sys.executable, "-c", "pass"], "fake")

    out = os.path.join(d, "hits.tsv")
    ok = P._run_besthit(Cfg(), t, os.path.join(d, "db"), q, out, "test")
    print("run ok            :", ok)

    rows = [l.split("\t") for l in open(out).read().splitlines() if l]
    names = sorted(r[0] for r in rows)
    expect = sorted("g%03d" % i for i in range(300))
    print("records           :", len(rows))
    print("every query once  :", names == expect)
    print("threads per shard :", sorted({r[1] for r in rows}))
    print("cache sidecar     :", os.path.exists(out + ".meta"))
    leftover = [f for f in os.listdir(d) if ".shard" in f or ".part" in f]
    print("shard leftovers   :", leftover)

    t0 = time.time()
    P._run_besthit(Cfg(), t, os.path.join(d, "db"), q, out, "test")
    print("2nd call (cached) : %.3fs" % (time.time() - t0))

    # a change to the inputs must invalidate the cache
    k1 = P._cache_key(Cfg(), t, q, "test")
    with open(q, "a") as o:
        o.write(">g999\nMKV\n")
    k2 = P._cache_key(Cfg(), t, q, "test")
    print("key changes when the query changes :", k1 != k2, "(must be True)")
    # ...and so must a change to the screening parameters
    class CfgHigh(Cfg):
        def i(self, sec, key, dv=0):
            if key == "max_target_seqs":
                return 999
            return Cfg.i(self, sec, key, dv)
    k3 = P._cache_key(CfgHigh(), t, q, "test")
    print("key changes when max_target_seqs changes :", k1 != k3, "(must be True)")

    shutil.rmtree(d)
    print("PASS" if names == expect and not leftover else "FAIL")


if __name__ == "__main__":
    main()
