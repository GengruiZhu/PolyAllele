#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PolyAllele v4 -- a graph-free, position-first pipeline for building standardised
allele matrices in polyploid genomes.

    Row   = one ancestral locus, defined by 1:1 collinear projection onto a
            diploid reference genome (or, failing that, by a reference-bracketed
            interval slot).
    Cell  = (locus, chromosome).  At most one gene is the representative; any
            further copies at the same locus are recorded explicitly as extras.

The grouping criterion is POSITION throughout.  Sequence similarity is used in
exactly two places, both of them confined to a window that position has already
authorised:
    * resolving which of several reference loci inside an interval a gene
      belongs to (stage `interp`, `--similarity-first`);
    * pulling a stray singleton back into a nearby rich row (stage `refine`).

Stages
    prep     reference split, target discovery, manifest, symlink farm
    project  pairwise JCVI collinear projection                 (tier 1)
    interp   order-preserving interpolation + synthetic slots   (tier 2)
    refine   position-constrained homology refinement           (tier 3)
    report   summary statistics
    all      prep -> project -> interp -> refine -> report

Usage
    python polyallele.py --config config.ini all
    python polyallele.py --config config.ini project      # resume a stage
"""

from __future__ import annotations

import argparse
import configparser
import csv
import glob
import hashlib
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict

VERSION = "4.4"
LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

# standard genetic code, for translating reference CDS when no protein file is given
CODON = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

NULLS = {"NA", "NAN", "", ".", "-"}


# --------------------------------------------------------------------------- #
#  logging
# --------------------------------------------------------------------------- #
class Log:
    def __init__(self, path=None):
        self.fh = open(path, "a", buffering=1) if path else None

    def __call__(self, *a):
        msg = " ".join(str(x) for x in a)
        stamp = time.strftime("[%Y-%m-%d %H:%M:%S]")
        line = f"{stamp} {msg}"
        print(line, flush=True)
        if self.fh:
            self.fh.write(line + "\n")

    def close(self):
        if self.fh:
            self.fh.close()


LOG = Log()


# --------------------------------------------------------------------------- #
#  configuration
# --------------------------------------------------------------------------- #
class Cfg:
    """Typed accessor over config.ini."""

    def __init__(self, path):
        self.path = os.path.abspath(path)
        cp = configparser.ConfigParser(inline_comment_prefixes=(";", "#"),
                                       interpolation=None)
        cp.optionxform = str.lower
        if not os.path.exists(path):
            raise SystemExit(f"config file not found: {path}")
        # utf-8-sig so a config written on Windows (which may carry a BOM) parses;
        # a BOM would otherwise make configparser reject line 1 as a missing
        # section header, with a message that does not mention the BOM at all.
        cp.read(path, encoding="utf-8-sig")
        self.cp = cp
        self.root = os.path.dirname(self.path)
        self._store = {}
        self.wd = ""
        self.od = ""
        # Every reader normalises gene IDs through this pattern, so the target
        # BED, the JCVI anchors and the GFF all agree on what a gene is called.
        set_id_strip(cp.get("input", "gene_id_strip", fallback=""))

    # inter-stage state (kept in memory during a run)
    def get(self, key, default=None):
        return self._store.get(key, default)

    def set(self, key, value):
        self._store[key] = value

    def raw(self, sec, key, default=""):
        try:
            return self.cp.get(sec, key).strip()
        except (configparser.NoSectionError, configparser.NoOptionError):
            return default

    def s(self, sec, key, default=""):
        v = self.raw(sec, key, default)
        return v if v else default

    def i(self, sec, key, default=0):
        v = self.raw(sec, key, "")
        try:
            return int(v)
        except ValueError:
            return default

    def f(self, sec, key, default=0.0):
        v = self.raw(sec, key, "")
        try:
            return float(v)
        except ValueError:
            return default

    def b(self, sec, key, default=False):
        v = self.raw(sec, key, "").lower()
        if v in ("1", "true", "yes", "on"):
            return True
        if v in ("0", "false", "no", "off"):
            return False
        return default

    def p(self, sec, key, default=""):
        """Path, resolved against the directory holding config.ini."""
        v = self.s(sec, key, default)
        if not v:
            return ""
        v = os.path.expanduser(os.path.expandvars(v))
        return v if os.path.isabs(v) else os.path.normpath(os.path.join(self.root, v))

    def list(self, sec, key, default=""):
        v = self.raw(sec, key, default)
        return [x.strip() for x in v.split(",") if x.strip()] if v else []


# --------------------------------------------------------------------------- #
#  small helpers
# --------------------------------------------------------------------------- #
def split_cell(v):
    v = str(v).strip()
    if v.upper() in NULLS:
        return []
    return [g.strip() for g in v.split(",") if g.strip() and g.strip().upper() not in NULLS]


# --------------------------------------------------------------------------- #
#  input readers
#
#  These are deliberately forgiving.  A gene ID has to match between the target
#  BED and the JCVI anchors or NOTHING is anchored, and every annotation pipeline
#  writes IDs a little differently -- a trailing `-P` / `-T` / `.1`, a quoted GFF
#  attribute, a CRLF from a Windows round trip.  A reader that insists on one
#  exact spelling turns all of that into a silent empty result.
# --------------------------------------------------------------------------- #
ID_STRIP = [None]          # compiled regex, or None; set from the config


def set_id_strip(pattern):
    """Compile `input.gene_id_strip` once."""
    pat = (pattern or "").strip()
    ID_STRIP[0] = re.compile(pat) if pat else None
    return ID_STRIP[0]


def norm_id(g):
    """Normalise a gene ID so the same gene written two ways still matches."""
    g = str(g).strip().strip('"').strip("'")
    rx = ID_STRIP[0]
    return rx.sub("", g) if rx else g


def read_fasta(path):
    name, buf = None, []
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line[0] in ";":            # legacy FASTA comment
                continue
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name = norm_id(line[1:].split()[0]) if line[1:].strip() else None
                buf = []
            elif name:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def read_bed(path):
    """Return [(start, gene), ...] sorted by position.

    Accepts tab- or space-separated BED, with or without track/browser/comment
    lines, and with or without a trailing column set.
    """
    out = []
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            p = line.split("\t") if "\t" in line else line.split()
            if len(p) >= 4:
                try:
                    out.append((int(p[1]), norm_id(p[3])))
                except ValueError:
                    continue
    out.sort()
    return out


def read_anchors(path):
    """JCVI .anchors: geneA <TAB> geneB <TAB> score."""
    out = []
    if not os.path.exists(path):
        return out
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if len(p) < 2:
                continue
            try:
                s = float(p[2])
            except (IndexError, ValueError):
                s = 0.0
            out.append((norm_id(p[0]), norm_id(p[1]), s))
    return out


def which(prog):
    return shutil.which(prog)


def translate(seq):
    seq = seq.upper().replace("U", "T")
    return "".join(CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def mkdir(path):
    os.makedirs(path, exist_ok=True)
    return path


def symlink(src, dst):
    if os.path.lexists(dst):
        return
    try:
        os.symlink(os.path.abspath(src), dst)
    except OSError:
        shutil.copy2(src, dst)


# --------------------------------------------------------------------------- #
#  STAGE 1 -- prep
# --------------------------------------------------------------------------- #
def parse_gff_genes(gff):
    """gene -> (chrom, start, end), from whatever the file actually contains.

    Annotation files disagree about almost everything:

      * some emit `gene` features, some only `mRNA` / `transcript` / `CDS`;
      * GFF3 uses `ID=`, GTF uses `gene_id=`, others use `Name=` or `locus_tag=`;
      * values may be quoted;
      * a Windows round trip leaves CRLF on every line.

    Reading only `gene` + `ID=` turns every other convention into zero genes, and
    the failure then surfaces much later as an empty run.  So: prefer `gene`
    features, count every attribute key that is known to carry an identifier, and
    fall back to the first gene-ish feature type when there are no `gene` lines.
    """
    ATTRS = ("ID", "gene_id", "Name", "gene", "locus_tag", "transcript_id",
             "Parent", "protein_id")
    FEAT_PREF = ("gene", "mRNA", "transcript", "CDS", "exon")
    by_type = defaultdict(dict)
    with open(gff, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 9:
                continue
            feat = p[2].strip()
            if feat not in FEAT_PREF:
                continue
            try:
                s, e = int(p[3]), int(p[4])
            except ValueError:
                continue
            attr = p[8]
            gid = None
            for key in ATTRS:
                m = re.search(rf"(?:^|;)\s*{key}\s*[= ]\s*"
                              r"(\"[^\"]+\"|[^;]+)", attr)
                if m:
                    gid = m.group(1)
                    break
            if not gid:
                continue
            gid = norm_id(gid)
            if gid and gid not in by_type[feat]:
                by_type[feat][gid] = (p[0], s, e)
    for feat in FEAT_PREF:                 # gene first, then the fallbacks
        if by_type.get(feat):
            return by_type[feat]
    return {}


def write_reference_frames(cfg, wd):
    """Split the reference into per-chromosome .cds / .bed / .pep."""
    ref = mkdir(os.path.join(wd, "reference"))
    gff = cfg.p("input", "reference_gff")
    cds = cfg.p("input", "reference_cds")
    pep = cfg.p("input", "reference_pep")
    if not gff or not cds:
        raise SystemExit("[prep] input.reference_gff and input.reference_cds are required")
    for f in (gff, cds):
        if not os.path.exists(f):
            raise SystemExit(f"[prep] missing: {f}")

    LOG(f"[prep] parsing reference GFF: {gff}")
    genes = parse_gff_genes(gff)
    bychrom = defaultdict(list)
    for g, (c, s, e) in genes.items():
        bychrom[c].append((s, e, g))
    LOG(f"[prep] reference: {len(genes):,} genes on {len(bychrom)} chromosomes")

    for c, lst in bychrom.items():
        lst.sort()
        with open(os.path.join(ref, c + ".bed"), "w") as o:
            for s, e, g in lst:
                o.write(f"{c}\t{s-1}\t{e}\t{g}\n")

    id2chrom = {g: v[0] for g, v in genes.items()}
    handles = {c: open(os.path.join(ref, c + ".cds"), "w") for c in bychrom}
    unmapped = 0
    with open(cds) as fh:
        cur = None
        for line in fh:
            if line.startswith(">"):
                gid = line[1:].split()[0]
                cur = handles.get(id2chrom.get(gid))
                if cur is None:
                    unmapped += 1
            if cur is not None:
                cur.write(line)
    for h in handles.values():
        h.close()
    LOG(f"[prep] reference CDS records not in GFF (skipped): {unmapped:,}")

    # proteins: use the supplied file when present, otherwise translate
    if pep and os.path.exists(pep):
        p2 = {}
        for name, seq in read_fasta(pep):
            p2.setdefault(name.split(".")[0], seq)
        n = 0
        with open(os.path.join(ref, "reference.pep"), "w") as o:
            for name, seq in read_fasta(pep):
                o.write(f">{name}\n{seq}\n")
                n += 1
        LOG(f"[prep] reference proteins copied from {pep}: {n:,}")
    else:
        n = 0
        with open(os.path.join(ref, "reference.pep"), "w") as o:
            for name, seq in read_fasta(cds):
                aa = translate(seq)
                if len(aa) >= 30:
                    o.write(f">{name}\n{aa}\n")
                    n += 1
        LOG(f"[prep] reference proteins translated from CDS: {n:,}")

    refchroms = sorted(bychrom)
    LOG(f"[prep] reference chromosomes: {refchroms}")
    return refchroms


def extract_target_from_scratch(cfg, wd, chroms):
    """Run gffread per chromosome to obtain .cds / .pep, and write .bed.

    Idempotent: a chromosome whose .cds already exists is skipped, so re-running
    the `prep` stage after an interruption is cheap.
    """
    gff = cfg.p("input", "target_gff")
    fa = cfg.p("input", "target_fasta")
    gffread = which("gffread")
    if not gffread:
        raise SystemExit("[prep] gffread not found in PATH "
                         "(needed to extract target CDS); install it or supply "
                         "input.target_prep_dir instead")
    tgt_dir = mkdir(os.path.join(wd, "target"))
    gff_dir = mkdir(os.path.join(tgt_dir, "gff"))
    fa_dir = mkdir(os.path.join(tgt_dir, "fa"))

    todo = [c for c in chroms
            if not (os.path.exists(os.path.join(tgt_dir, c + ".cds"))
                    and os.path.getsize(os.path.join(tgt_dir, c + ".cds")) > 0)]
    if not todo:
        LOG(f"[prep] target sequences already extracted ({len(chroms)} chromosomes)")
        return tgt_dir

    LOG(f"[prep] extracting target sequences with gffread "
        f"({len(todo)} of {len(chroms)} chromosomes)")
    genes = parse_gff_genes(gff)
    per = defaultdict(list)
    for g, (c, s, e) in genes.items():
        if c in chroms:
            per[c].append((s, e, g))

    # index the genome FASTA by contig, then slice per chromosome
    want = set(todo)
    cur = None
    with open(fa) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    cur.close()
                    cur = None
                cid = line[1:].split()[0]
                if cid in want:
                    cur = open(os.path.join(fa_dir, cid + ".fa"), "w")
                    cur.write(line)
            elif cur:
                cur.write(line)
    if cur:
        cur.close()

    for c in sorted(todo):
        lst = sorted(per.get(c, []))
        with open(os.path.join(tgt_dir, c + ".bed"), "w") as o:
            for s, e, g in lst:
                o.write(f"{c}\t{s-1}\t{e}\t{g}\n")
        with open(os.path.join(gff_dir, c + ".gff"), "w") as o:
            o.write("##gff-version 3\n")
            with open(gff) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    if line.split("\t", 1)[0] == c:
                        o.write(line)
        subprocess.run([gffread, os.path.join(gff_dir, c + ".gff"),
                        "-g", os.path.join(fa_dir, c + ".fa"),
                        "-x", os.path.join(tgt_dir, c + ".cds"),
                        "-y", os.path.join(tgt_dir, c + ".pep")],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    LOG(f"[prep] gffread finished for {len(todo)} chromosomes")
    return tgt_dir


def resolve_groups(cfg):
    """Return [(group_name, ref_chrom, [chromosome, ...]), ...]."""
    mode = cfg.s("chromosomes", "mode", "group_file").lower()
    ref_prefix = cfg.s("reference", "prefix", "Chr")
    refchroms = cfg.get("_refchroms", [])

    def auto_ref(gname):
        m = re.search(r"(\d+)", gname)
        if not m:
            return ""
        num = int(m.group(1))
        for rc in refchroms:
            mm = re.search(r"(\d+)$", rc)
            if mm and int(mm.group(1)) == num:
                return rc
        return ""

    if mode == "group_file":
        gf = cfg.p("chromosomes", "group_file")
        if not gf or not os.path.exists(gf):
            raise SystemExit(f"[prep] chromosomes.group_file not found: {gf}")
        groups = []
        with open(gf) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                p = [x.strip() for x in line.rstrip("\n").split("\t") if x.strip()]
                if len(p) >= 2:
                    groups.append([p[0], p[1:]])
    else:
        prefix = cfg.s("chromosomes", "prefix", "Chr")
        g0 = cfg.i("chromosomes", "group_start", 1)
        g1 = cfg.i("chromosomes", "group_end", 10)
        subs = cfg.list("chromosomes", "subgenomes", "A,B")
        allchroms = cfg.get("_targetchroms", [])
        groups = []
        for g in range(g0, g1 + 1):
            chroms = [f"{prefix}{g}{s}" for s in subs]
            chroms = [c for c in chroms if c in allchroms]
            if chroms:
                groups.append([f"Group_{g:02d}", chroms])

    mapping = {}
    for pair in cfg.list("reference", "group_ref_map"):
        if ":" in pair:
            k, v = pair.split(":", 1)
            mapping[k.strip()] = v.strip()

    out = []
    for gname, chroms in groups:
        ref = mapping.get(gname) or mapping.get(re.sub(r"\D", "", gname)) or auto_ref(gname)
        if not ref:
            raise SystemExit(f"[prep] cannot determine the reference chromosome for {gname}; "
                             f"set reference.group_ref_map")
        if ref not in refchroms:
            raise SystemExit(f"[prep] {gname} -> {ref} is not a reference chromosome "
                             f"({refchroms})")
        out.append((gname, ref, chroms))
    return out


def stage_prep(cfg):
    wd = cfg.wd
    tgt = mkdir(os.path.join(wd, "target"))
    farm = mkdir(os.path.join(wd, "farm"))

    refchroms = write_reference_frames(cfg, wd)
    cfg.set("_refchroms", refchroms)

    prep_dir = cfg.p("input", "target_prep_dir")
    if prep_dir:
        if not os.path.isdir(prep_dir):
            raise SystemExit(f"[prep] input.target_prep_dir not found: {prep_dir}")
        all_cds = sorted(glob.glob(os.path.join(prep_dir, "**", "*.cds"), recursive=True))
        cds_files, seen = [], set()
        for f in all_cds:
            c = os.path.basename(f)[:-4]
            if c in seen:
                LOG(f"[prep]   duplicate CDS for {c}, ignoring {f}")
                continue
            seen.add(c)
            cds_files.append(f)
        chroms = [os.path.basename(f)[:-4] for f in cds_files]
        LOG(f"[prep] target_prep_dir mode: {len(chroms)} chromosomes found")
        # .cds / .bed / .pep are often in sibling directories, so look them up
        # by basename anywhere under target_prep_dir
        pool = {}
        for ext in (".bed", ".pep"):
            for f in glob.glob(os.path.join(prep_dir, "**", "*" + ext), recursive=True):
                pool.setdefault(os.path.basename(f), f)
        derived = None
        for f in cds_files:
            c = os.path.basename(f)[:-4]
            symlink(f, os.path.join(tgt, c + ".cds"))
            for ext in (".bed", ".pep"):
                src = pool.get(c + ext)
                if src:
                    symlink(src, os.path.join(tgt, c + ext))
            if not os.path.exists(os.path.join(tgt, c + ".bed")):
                if derived is None and cfg.p("input", "target_gff"):
                    LOG("[prep] no BED files found; deriving them from input.target_gff")
                    derived = parse_gff_genes(cfg.p("input", "target_gff"))
                if derived:
                    with open(os.path.join(tgt, c + ".bed"), "w") as o:
                        for s, e, g in sorted(v for v in
                                              ((s, e, g) for g, (cc, s, e) in derived.items()
                                               if cc == c)):
                            o.write(f"{c}\t{s-1}\t{e}\t{g}\n")
                else:
                    LOG(f"[prep] WARNING: no BED for {c} "
                        f"(set input.target_gff so it can be derived)")
    else:
        if not cfg.p("input", "target_gff") or not cfg.p("input", "target_fasta"):
            raise SystemExit("[prep] give either input.target_prep_dir, or both "
                             "input.target_gff and input.target_fasta")
        genes = parse_gff_genes(cfg.p("input", "target_gff"))
        chroms = sorted({v[0] for v in genes.values()})
        tgt_dir = extract_target_from_scratch(cfg, wd, chroms)
        for c in chroms:
            for ext in (".cds", ".bed", ".pep"):
                src = os.path.join(tgt_dir, c + ext)
                if os.path.exists(src):
                    symlink(src, os.path.join(tgt, c + ext))
    cfg.set("_targetchroms", chroms)

    # Fail here, with something actionable, rather than 40 minutes later with an
    # empty matrix.  Every way this can go wrong looks the same from downstream.
    if not chroms:
        raise SystemExit(
            "[prep] no target chromosomes were found.\n"
            "       * with input.target_prep_dir: the directory must contain\n"
            "         <chromosome>.cds files (searched recursively); or\n"
            "       * with input.target_gff: the file must hold gene / mRNA /\n"
            "         transcript / CDS records in GFF3 or GTF form.\n"
            "       If the records are there but nothing was read, the identifier\n"
            "       attribute or the ID spelling is the likely cause -- see\n"
            "       input.gene_id_strip in the [input] section.")

    groups = resolve_groups(cfg)
    cfg.set("_groups", groups)

    man = os.path.join(wd, "manifest.tsv")
    with open(man, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["Group", "RefChrom", "Chromosome", "CDS", "BED", "PEP"])
        for gname, ref, gc in groups:
            for c in gc:
                cds = os.path.join(tgt, c + ".cds")
                if not os.path.exists(cds):
                    raise SystemExit(f"[prep] no CDS for target chromosome {c} ({cds})")
                w.writerow([gname, ref, c, cds,
                            os.path.join(tgt, c + ".bed"),
                            os.path.join(tgt, c + ".pep")])

    # symlink farm for the JCVI working directories
    for d in (os.path.join(wd, "jcvi"), os.path.join(wd, "jcvi_audit")):
        mkdir(d)
        for f in glob.glob(os.path.join(tgt, "*")):
            symlink(f, os.path.join(d, os.path.basename(f)))
        for f in glob.glob(os.path.join(wd, "reference", "*.cds")) + \
                 glob.glob(os.path.join(wd, "reference", "*.bed")):
            symlink(f, os.path.join(d, os.path.basename(f)))

    LOG("")
    LOG(f"[prep] {len(groups)} homology groups, {len(chroms)} target chromosomes")
    for gname, ref, gc in groups:
        LOG(f"[prep]   {gname} -> reference {ref}: {len(gc):2d} chromosomes "
            f"{gc[0]}..{gc[-1]}")
    LOG(f"[prep] manifest: {man}")


# --------------------------------------------------------------------------- #
#  STAGE 2 -- projection (JCVI)
# --------------------------------------------------------------------------- #
def run_pair(cfg, d, a, b, cscore):
    out = os.path.join(d, f"{a}.{b}.anchors")
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return "skip"
    cpus = []
    try:
        h = subprocess.run([sys.executable, "-m", "jcvi.compara.catalog", "ortholog", "--help"],
                           capture_output=True, text=True)
        if "--cpus" in (h.stdout + h.stderr):
            cpus = ["--cpus", str(cfg.i("projection", "threads", 8))]
    except Exception:
        pass
    log = os.path.join(d, f"{a}.{b}.jcvi.log")
    with open(log, "w") as lf:
        rc = subprocess.run([sys.executable, "-m", "jcvi.compara.catalog", "ortholog",
                             "--no_strip_names", a, b, f"--cscore={cscore}", "--no_dotplot"]
                            + cpus, cwd=d, stdout=lf, stderr=lf).returncode
    return "ok" if rc == 0 else f"fail:{rc}"


def stage_project(cfg):
    wd = cfg.wd
    cscore = cfg.f("projection", "cscore", 0.70)
    groups = cfg.get("_groups")
    d = os.path.join(wd, "jcvi")
    pairs = [(c, ref) for _g, ref, gc in groups for c in gc]

    # Optional: reuse anchors that were computed earlier (e.g. by a previous run
    # or by another pipeline), instead of recomputing the pairwise alignments.
    reuse = cfg.p("input", "anchors_dir")
    pool = {}
    if reuse:
        if not os.path.isdir(reuse):
            raise SystemExit(f"[project] input.anchors_dir not found: {reuse}")
        for f in glob.glob(os.path.join(reuse, "**", "*.anchors"), recursive=True):
            pool.setdefault(os.path.basename(f), f)
        LOG(f"[project] reusing precomputed anchors from {reuse} "
            f"({len(pool):,} .anchors files available)")

    LOG(f"[project] {len(pairs)} target-vs-reference pairs at cscore={cscore}")
    t0 = time.time()
    tally = Counter()
    for i, (a, b) in enumerate(pairs, 1):
        src = pool.get(f"{a}.{b}.anchors") if pool else None
        if src:
            symlink(src, os.path.join(d, f"{a}.{b}.anchors"))
            tally["reused"] += 1
            continue
        r = run_pair(cfg, d, a, b, cscore)
        tally[r.split(":")[0]] += 1
        if r.startswith("fail"):
            LOG(f"[project]   FAIL {a} vs {b} ({r}) -> {a}.{b}.jcvi.log")
        elif i % 20 == 0:
            LOG(f"[project]   {i}/{len(pairs)} done, {int(time.time()-t0)}s")
    LOG(f"[project] finished in {int(time.time()-t0)}s: "
        + ", ".join(f"{k}={v}" for k, v in sorted(tally.items())))

    if cfg.b("projection", "audit", False):
        da = os.path.join(wd, "jcvi_audit")
        refchroms = cfg.get("_refchroms")
        LOG(f"[project] audit: {len(pairs)}x{len(refchroms)} pairs "
            f"(identity criterion, for comparison only)")
        n = 0
        for a, b in pairs:
            for rc in refchroms:
                run_pair(cfg, da, a, rc, cscore)
                n += 1
                if n % 100 == 0:
                    LOG(f"[project]   audit {n} pairs")
        LOG("[project] audit finished")


# --------------------------------------------------------------------------- #
#  STAGE 3 -- interpolation (tier 1 + tier 2)
# --------------------------------------------------------------------------- #
def load_manifest(wd):
    man = []
    with open(os.path.join(wd, "manifest.tsv")) as fh:
        for x in csv.DictReader(fh, delimiter="\t"):
            man.append((x["Group"], x["RefChrom"], x["Chromosome"], x["CDS"], x["PEP"]))
    return man


def load_reference_order(wd):
    order, lists = {}, {}
    for bed in sorted(glob.glob(os.path.join(wd, "reference", "*.bed"))):
        c = os.path.basename(bed)[:-4]
        genes = [g for _s, g in read_bed(bed)]
        lists[c] = genes
        for i, g in enumerate(genes):
            order[g] = (c, i)
    return order, lists


def interpolate(chrom, genes, tier1, eru_order, eru_list, max_gap,
                max_same_locus=50, keep_collisions=True, tandem_window=5,
                distribute_surplus=False):
    """gene -> (locus, source, interval, synthetic, bounds).

    Two situations get their own comment, because both were once handled wrongly
    and manufactured spurious singleton loci:

    1. Both flanking anchors sit on the SAME reference locus.  The genes between
       them are bracketed evidence that they are additional copies AT that locus.

    2. Two genes on the same chromosome project to the same reference locus.
       A genuine tandem array looks exactly like this -- and so does a projection
       artefact, because intervals in reference-locus space can overlap when the
       projection is non-monotonic.  The two are told apart by GENE ORDER: a real
       array puts neighbours together, an artefact puts genes that are far apart
       along the chromosome into one cell.

       Only neighbours within `tandem_window` gene models are therefore kept in
       the cell as extra copies.  Anything further away gets its own slot, which
       is what the pipeline did before v4.2 -- and which the third-party
       comparison showed to be the safer choice: keeping every collision doubled
       the provably-bad merge rate (6.54% -> 15.12%) and dropped agreement with an
       independent method from 81.4% to 66.5%.
    """
    res = {}
    occupied = {}                 # locus -> [(gene, gene_order_index), ...]
    for i, g in enumerate(genes):
        if g in tier1 and tier1[g][0] in eru_order:
            res[g] = (tier1[g][0], "anchor", "", False, None)
            occupied.setdefault(tier1[g][0], []).append((g, i))
    synth_ctr = defaultdict(int)

    def slot(iv, tag):
        k = synth_ctr[(iv, tag)]
        synth_ctr[(iv, tag)] = k + 1
        return f"{iv}|{tag}{k}"

    pos_idx = [i for i, g in enumerate(genes) if g in res]
    stats = Counter()
    unplaced = {}            # gene -> (interval, bounds, reason)
    for a, b in zip(pos_idx, pos_idx[1:]):
        if b - a <= 1:
            continue
        o_u, o_d = eru_order.get(res[genes[a]][0]), eru_order.get(res[genes[b]][0])
        if not o_u or not o_d or o_u[0] != o_d[0]:
            for idx in range(a + 1, b):
                unplaced[genes[idx]] = ("", None,
                                        "flanks_on_different_reference_chromosomes")
            stats["flanks_on_different_reference_chromosomes"] += b - a - 1
            continue
        i, j = o_u[1], o_d[1]
        mids = list(range(a + 1, b))

        # ---- both anchors on the SAME reference locus -----------------------
        if i == j:
            L = eru_list[o_u[0]][i]
            iv = f"{o_u[0]}|same|{i}"
            for k, idx in enumerate(mids):
                g = genes[idx]
                if k < max_same_locus:
                    res[g] = (L, "copy_at_anchored_locus", iv, False, (i, i))
                    occupied.setdefault(L, []).append((g, idx))
                    stats["copy_at_anchored_locus"] += 1
                else:
                    res[g] = (slot(iv, "c"), "copy_at_anchored_locus_over_cap",
                              iv, True, (i, i))
                    stats["copy_at_anchored_locus_over_cap"] += 1
            continue

        inv = j < i
        lo, hi = (j, i) if inv else (i, j)
        g_gap = hi - lo - 1
        m = len(mids)
        if g_gap > max_gap:
            iv = f"{o_u[0]}|wide|{lo}-{hi}"
            for idx in mids:
                unplaced[genes[idx]] = (iv, (lo, hi),
                                        "interval_wider_than_max_interval")
            stats["interval_wider_than_max_interval"] += m
            continue
        sub = eru_list[o_u[0]][lo + 1:hi]
        if inv:
            sub = sub[::-1]
        iv = f"{o_u[0]}|iv|{lo}-{hi}" + ("inv" if inv else "")

        def place(idx, cand, src):
            g = genes[idx]
            # a slot id carries |iv| or |same|; a real locus is a reference gene
            syn = ("|iv|" in cand) or ("|same|" in cand)
            cell = occupied.get(cand)
            if cell is None:
                occupied[cand] = [(g, idx)]
                res[g] = (cand, src, iv, syn, (lo, hi))
                stats[src] += 1
                return
            # Collisions are only absorbed when the gene is a NEIGHBOUR of a gene
            # already in that cell: that is what a real tandem array looks like.
            # A collision between genes far apart along the chromosome is a
            # projection artefact (overlapping intervals), and merging them would
            # put two different loci in one cell.
            near = any(abs(idx - j) <= tandem_window for _gj, j in cell)
            if keep_collisions and near and len(cell) < max_same_locus:
                cell.append((g, idx))
                res[g] = (cand, src.replace("_copy", "") + "_copy", iv, False, (lo, hi))
                stats["collision_kept_as_tandem_copy"] += 1
            else:
                L = slot(iv, "s")
                occupied[L] = [(g, idx)]
                res[g] = (L, src + "_pushed", iv, True, (lo, hi))
                stats["collision_pushed_to_synthetic"] += 1

        if g_gap == 0:
            # the two anchors are adjacent reference loci: no room at all
            for idx in mids:
                place(idx, slot(iv, "n"), "interp_synthetic")
        elif m == g_gap:
            for k, idx in enumerate(mids):
                place(idx, sub[k], "interp_exact" + ("_inv" if inv else ""))
        elif m < g_gap:
            for k, idx in enumerate(mids):
                t = (g_gap // 2) if m == 1 else int(round(k * (g_gap - 1) / (m - 1)))
                place(idx, sub[t], "interp_proportional" + ("_inv" if inv else ""))
        else:
            # More target genes than reference loci in this interval.
            #
            # `distribute_surplus = true` spreads them over the interval's loci,
            # calling the surplus extra copies.  That sounds helpful -- it avoids
            # manufacturing singleton rows -- but attributing the third-party
            # disagreements back to their `Source` showed this branch produced
            # about HALF of all provably-bad merges, more than any other stage.
            # The reason is simple: when the reference offers four loci and the
            # target has ten genes, calling six of them "copies at those loci" is
            # an inference, not evidence, and the position of those six is exactly
            # what is unknown.
            #
            # The default is therefore the conservative reading: the first g_gap
            # genes take the loci in order and the rest get their own slot.
            for k, idx in enumerate(mids):
                if distribute_surplus:
                    t = min(int(k * g_gap / m), g_gap - 1)
                    place(idx, sub[t], "interp_distributed" + ("_inv" if inv else ""))
                elif k < g_gap:
                    place(idx, sub[k], "interp_exact" + ("_inv" if inv else ""))
                else:
                    place(idx, slot(iv, "x"), "interp_synthetic")
    head = pos_idx[0] if pos_idx else len(genes)
    for idx in range(head):
        unplaced[genes[idx]] = ("", None, "before_the_first_anchor")
    stats["before_the_first_anchor"] += head
    tail_start = (pos_idx[-1] + 1) if pos_idx else 0
    for idx in range(tail_start, len(genes)):
        unplaced[genes[idx]] = ("", None, "after_the_last_anchor")
    stats["after_the_last_anchor"] += len(genes) - tail_start
    return res, stats, unplaced


def stage_interp(cfg):
    wd = cfg.wd
    man = load_manifest(wd)
    eru_order, eru_list = load_reference_order(wd)
    global _REFGENES
    _REFGENES = set(eru_order)
    cfg.set("_eru_order", eru_order)
    cfg.set("_eru_lists", eru_list)
    LOG(f"[interp] reference loci: {len(eru_order):,}")

    tier1 = {}
    for _g, ref, chrom, _c, _p in man:
        tg = defaultdict(list)
        for a, b, s in read_anchors(os.path.join(wd, "jcvi", f"{chrom}.{ref}.anchors")):
            if a in _REFGENES and b not in _REFGENES:
                tg[b].append((a, s))
            elif b in _REFGENES and a not in _REFGENES:
                tg[a].append((b, s))
        for g, lst in tg.items():
            tier1[g] = max(lst, key=lambda x: x[1])
    LOG(f"[interp] tier-1 anchored genes: {len(tier1):,}")

    max_gap = cfg.i("interpolation", "max_interval", 50)
    max_same = cfg.i("interpolation", "max_copies_at_anchored_locus", 50)
    keep_coll = cfg.b("interpolation", "keep_same_chromosome_collisions", True)
    tandem_win = cfg.i("interpolation", "tandem_gene_window", 5)
    distribute = cfg.b("interpolation", "distribute_surplus", False)
    gene_loc, gene_meta = {}, {}
    unplaced = {}
    allstats = Counter()
    totals = Counter()
    for gname, ref, chrom, cds, pep in man:
        genes = [g for _s, g in read_bed(os.path.join(wd, "target", chrom + ".bed"))]
        totals["genes"] += len(genes)
        res, st, unp = interpolate(chrom, genes, tier1, eru_order, eru_list,
                                   max_gap, max_same, keep_coll, tandem_win, distribute)
        allstats.update(st)
        for g, (L, src, iv, syn, bnd) in res.items():
            gene_loc[g] = L
            gene_meta[g] = (src, iv, gname, chrom, syn, bnd)
        for g, (iv, bnd, why) in unp.items():
            unplaced[g] = (iv, bnd, why, gname, chrom)
        LOG(f"[interp] {chrom:12s} genes={len(genes):7,} placed={len(res):7,} "
            f"({100.0*len(res)/max(len(genes),1):5.2f}%)")

    LOG("")
    LOG("[interp] assignment cases:")
    for k, v in sorted(allstats.items()):
        LOG(f"[interp]   {k:48s} {v:>9,}")

    cfg.set("_gene_loc", gene_loc)
    cfg.set("_gene_meta", gene_meta)
    cfg.set("_tier1", tier1)
    cfg.set("_unplaced", unplaced)

    # optional: BLASTP disambiguation of proportional assignments inside intervals
    if cfg.b("refinement", "similarity_first", True):
        _similarity_first(cfg, man)

def _similarity_first(cfg, man):
    """For genes placed by proportional spacing, use BLASTP inside the interval."""
    gene_loc = cfg.get("_gene_loc")
    gene_meta = cfg.get("_gene_meta")
    cands = [g for g, m in gene_meta.items() if m[0].startswith("interp_proportional")]
    if not cands:
        return
    if not (which("diamond") or (which("blastp") and which("makeblastdb"))):
        LOG("[interp] similarity pass skipped: no aligner in PATH")
        return
    wd = cfg.wd
    tmp = mkdir(os.path.join(wd, "tmp"))
    LOG(f"[interp] similarity pass for {len(cands):,} proportional assignments")

    pep_of = {}
    for _g, _ref, chrom, _c, pep in man:
        pep_of[chrom] = pep
    sfa = os.path.join(wd, "reference", "reference.pep")
    qfa = os.path.join(tmp, "query.fa")
    n = _write_pep(pep_of, cands, qfa)
    db = os.path.join(tmp, "refpep")
    bo = os.path.join(tmp, "similarity_first.tsv")
    try:
        if not _run_besthit(cfg, sfa, db, qfa, bo, "interp"):
            return
    except Exception as e:
        LOG(f"[interp] similarity pass failed ({e}); keeping proportional placements")
        return

    # Eru gene -> allowed interval set per query
    allowed = defaultdict(set)
    for g in cands:
        bnd = gene_meta[g][5]
        iv = gene_meta[g][1]
        if bnd:
            allowed[g] = (iv.split("|iv|")[0], bnd[0], bnd[1])
    best = {}
    with open(bo) as fh:
        for line in fh:
            p = line.split("\t")
            if len(p) < 5:
                continue
            q, s, pid, ln, bs = p[0], p[1], float(p[2]), int(p[3]), float(p[4])
            a = allowed.get(q)
            if not a:
                continue
            gpos = cfg.get("_eru_order").get(s)
            if not gpos or gpos[0] != a[0] or not (a[1] <= gpos[1] <= a[2]):
                continue
            if q not in best or bs > best[q][1]:
                best[q] = (s, bs)
    nfix = 0
    for q, (s, _bs) in best.items():
        old = gene_loc[q]
        if old == s:
            continue
        gene_loc[q] = s
        src, iv, gn, ch, _sy, bnd = gene_meta[q]
        gene_meta[q] = (src + "_similarity", iv, gn, ch, False, bnd)
        nfix += 1
    LOG(f"[interp] similarity pass reassigned {nfix:,} of {len(cands):,} genes")


# --------------------------------------------------------------------------- #
#  STAGE 3.5 -- locus repair
# --------------------------------------------------------------------------- #
def _write_pep(pep_of, genes, path):
    """Write a FASTA holding exactly `genes`, in the order the files offer."""
    want = set(genes)
    n = 0
    with open(path, "w") as o:
        for chrom, pf in sorted(pep_of.items()):
            if not want or not os.path.exists(pf):
                continue
            for name, seq in read_fasta(pf):
                if name in want:
                    o.write(f">{name}\n{seq}\n")
                    want.discard(name)
                    n += 1
                    if not want:
                        break
    return n


def _search_cmd(cfg, used, query_fa, db_prefix, out_tsv, threads, label=""):
    """The search half of the best-hit screen, for one query file.

    `max_target_seqs` is resolved per stage, because the three screens do not
    want the same thing.  `interp` and `refine` pick the best hit inside a small
    candidate window and were validated at 20; `repair` searches the whole placed
    set for an allelic partner and needs a wider net (50).  A per-stage key wins
    when present:
        max_target_seqs_<label>   e.g. max_target_seqs_repair
        max_target_seqs           the fallback
    """
    ev = cfg.s("alignment", "evalue", "1e-5")
    kmax = cfg.i("alignment", f"max_target_seqs_{label}", 0) if label else 0
    if not kmax:
        kmax = cfg.i("alignment", "max_target_seqs", 20)
    kmax = str(kmax)
    fields = "qseqid sseqid pident length bitscore"
    if used == "diamond":
        return (["diamond", "blastp", "-q", query_fa, "-d", db_prefix,
                 "-o", out_tsv, "-f", "6"] + fields.split()
                + ["-e", ev, "-k", kmax, "--threads", str(threads), "--sensitive"])
    return ["blastp", "-query", query_fa, "-db", db_prefix,
            "-outfmt", "6 " + fields, "-evalue", ev, "-max_target_seqs", kmax,
            "-num_threads", str(threads), "-out", out_tsv]


def _besthit_cmds(cfg, target_fa, db_prefix, query_fa, out_tsv):
    """Build the database for the best-hit screen and pick the aligner.

    blastp by default.  DIAMOND is offered as `runtime.aligner = diamond` because
    it is about two orders of magnitude faster, but it is NOT a drop-in: measured
    on 2,000 queries against the same candidate set it returned the same best
    partner as blastp for only 58.65% of them, and the sets passing the pipeline's
    thresholds differed by ~14% (Jaccard 86.16% at id>=70/len>=60, 85.72% at
    id>=90/len>=80).  In a polyploid many homoeologs are near-identical, so top
    hits are near-ties and the two aligners break them differently.  See config.ini.
    """
    choice = cfg.s("runtime", "aligner", "blastp").lower()
    if choice == "diamond" and which("diamond"):
        return (["diamond", "makedb", "--in", target_fa, "-d", db_prefix], "diamond")
    if not which("makeblastdb") or not which("blastp"):
        return None, "none"
    return (["makeblastdb", "-in", target_fa, "-dbtype", "prot", "-out", db_prefix],
            "blastp")


def _cache_key(cfg, target_fa, query_fa, label):
    """Everything that decides what the screen would produce.

    Sizes alone are not enough: two runs can screen the same files with different
    `evalue`, `max_target_seqs` or aligner and get different hits.  Including a
    digest of the search parameters is what keeps the cache honest -- without it
    a run configured for k=20 silently reused a cached k=50 result.
    """
    parts = [str(os.path.getsize(query_fa)), str(os.path.getsize(target_fa))]
    for sec, key in (("alignment", "evalue"),
                     ("alignment", f"max_target_seqs_{label}"),
                     ("alignment", "max_target_seqs"),
                     ("runtime", "aligner")):
        parts.append(f"{key}={cfg.s(sec, key, '')}")
    return hashlib.sha1("|".join(parts).encode()).hexdigest()


def _cache_valid(out_tsv, key):
    meta = out_tsv + ".meta"
    if not (os.path.exists(out_tsv) and os.path.getsize(out_tsv) > 0
            and os.path.exists(meta)):
        return False
    try:
        with open(meta) as fh:
            return fh.read().strip() == key
    except OSError:
        return False


def _split_fasta(path, n, outdir, tag):
    """Split a FASTA into n roughly equal sequential chunks.

    The shard is chosen from the byte offset at which a record STARTS, so `seen`
    must accumulate EVERY line -- headers and sequence alike.  Counting only the
    headers makes `seen` grow an order of magnitude slower than the file, which
    silently puts every record in the first shard.
    """
    if n <= 1:
        return [path]
    paths = [os.path.join(outdir, f".{tag}.shard{i}.pep") for i in range(n)]
    handles = [open(p, "w") for p in paths]
    total = max(os.path.getsize(path), 1)
    idx, seen = 0, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                idx = min(int(seen * n / total), n - 1)
            seen += len(line)
            handles[idx].write(line)
    for h in handles:
        h.close()
    return paths


def _run_besthit(cfg, target_fa, db_prefix, query_fa, out_tsv, label):
    build, used = _besthit_cmds(cfg, target_fa, db_prefix, query_fa, out_tsv)
    if used == "none":
        LOG(f"[{label}] skipped: neither diamond nor blastp/makeblastdb in PATH")
        return False
    ckey = _cache_key(cfg, target_fa, query_fa, label)
    if not cfg.b("alignment", "force", False) and _cache_valid(out_tsv, ckey):
        LOG(f"[{label}] reusing cached {os.path.basename(out_tsv)} "
            f"({os.path.getsize(out_tsv)/1e6:.1f} MB) -- delete it or set "
            f"alignment.force = true to redo")
        return True

    # Sharding.  blastp scales poorly with thread count, so on a node with many
    # cores it is faster to run several query shards side by side than one wide
    # process.  The aligner, its parameters, and therefore the resulting hits, are
    # unchanged; only the wall time is.
    threads = cfg.i("projection", "threads", 8)
    shards = max(1, min(cfg.i("alignment", "shards", 1), max(threads // 2, 1)))
    tmpdir = os.path.dirname(out_tsv) or "."
    qfiles = _split_fasta(query_fa, shards, tmpdir, label)
    per = max(1, threads // len(qfiles))
    errlog = os.path.join(tmpdir, f"{label}.aligner.err")

    LOG(f"[{label}] screening with {used}"
        + (f", {len(qfiles)} shards x {per} threads" if len(qfiles) > 1
           else f", {per} threads"))
    cmd0 = _search_cmd(cfg, used, qfiles[0], db_prefix, out_tsv, per, label)
    LOG(f"[{label}]   " + " ".join(str(x) for x in cmd0[:8])
        + (" ..." if len(cmd0) > 8 else ""))
    t0 = time.time()
    errf = open(errlog, "w")
    try:
        subprocess.run(build, check=True, stdout=subprocess.DEVNULL, stderr=errf)
        parts = [out_tsv + f".part{i}" for i in range(len(qfiles))]
        procs = [subprocess.Popen(
            _search_cmd(cfg, used, q, db_prefix, o, per, label),
            stdout=subprocess.DEVNULL, stderr=errf) for q, o in zip(qfiles, parts)]
        bad = [p.returncode for p in procs if p.wait() != 0]
        if bad:
            raise subprocess.CalledProcessError(bad[0], used)
    except subprocess.CalledProcessError as e:
        errf.close()
        LOG(f"[{label}] {used} FAILED (exit {e.returncode}); see {errlog}")
        try:
            with open(errlog) as ef:
                for line in ef.readlines()[-8:]:
                    LOG(f"[{label}]   | {line.rstrip()}")
        except OSError:
            pass
        raise
    errf.close()
    with open(out_tsv, "w") as o:
        for pf in parts:
            if os.path.exists(pf):
                with open(pf) as fh:
                    shutil.copyfileobj(fh, o)
                os.remove(pf)
    for q in qfiles:
        if q != query_fa and os.path.exists(q):
            os.remove(q)

    try:
        sz = os.path.getsize(out_tsv)
    except OSError:
        sz = 0
    if sz > 0:
        with open(out_tsv + ".meta", "w") as fh:
            fh.write(ckey + "\n")
    LOG(f"[{label}] {used} finished in {int(time.time()-t0)}s, "
        f"{sz/1e6:.1f} MB of hits")
    return True


def _read_hits(path):
    hits = defaultdict(list)
    if not os.path.exists(path):
        return hits
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                try:
                    hits[p[0]].append((p[1], float(p[2]), int(p[3]), float(p[4])))
                except ValueError:
                    continue
    return hits


def stage_repair(cfg):
    """Similarity-guided locus repair.

    A gene the projection could not put at a reference locus -- it sits in a
    synthetic slot, or fell into an interval too wide to interpolate -- may still
    have an allele that WAS placed: its homoeolog on another chromosome of the
    same homology group.

    The move is allowed only when ALL of these hold:

      * the partner sits at a REAL reference locus L (not a synthetic slot);
      * same homology group but a DIFFERENT chromosome.  Genes on homoeologous
        chromosomes are allelic by default; a same-chromosome partner would be a
        paralog, and that case is refused outright;
      * BLASTP identity >= min_identity over >= min_alignment amino acids;
      * L lies inside the interval the gene was bracketed by, allowing
        `repair_window` reference loci of slack.

    Position authorises; similarity only chooses among what position allows.
    """
    if not cfg.b("repair", "enabled", True):
        LOG("[repair] disabled in config")
        return
    wd = cfg.wd
    gene_loc = cfg.get("_gene_loc")
    gene_meta = cfg.get("_gene_meta")
    unplaced = cfg.get("_unplaced") or {}
    eru_order = cfg.get("_eru_order") or {}
    man = load_manifest(wd)

    MINID = cfg.f("repair", "min_identity", 90)
    MINLEN = cfg.i("repair", "min_alignment", 80)
    WIN = cfg.i("repair", "window", 1)

    # ---- who may be repaired ------------------------------------------------
    # A candidate is either a gene sitting in a synthetic slot, or one that
    # entered no row at all but WAS bracketed by two anchors.  The two live in
    # different dicts, and a gene from the second one has no gene_meta entry yet,
    # so carry everything the move needs in the tuple.
    query = {}                      # gene -> (group, chrom, (lo,hi), interval)
    for g, m in gene_meta.items():
        src, iv, gn, ch, syn, bnd = m
        if syn or src.endswith("_pushed"):
            query[g] = (gn, ch, bnd, iv)
    for g, (iv, bnd, why, gn, ch) in unplaced.items():
        if bnd is not None:                       # only the bracketed ones
            query[g] = (gn, ch, bnd, iv)
    if not query:
        LOG("[repair] nothing to repair")
        return
    LOG(f"[repair] candidates            : {len(query):,} "
        f"(synthetic slots + wide intervals)")

    # ---- who may serve as the target ---------------------------------------
    target = {}                     # gene -> (group, chrom, locus, locus index)
    for g, m in gene_meta.items():
        src, iv, gn, ch, syn, bnd = m
        if syn:
            continue
        L = gene_loc.get(g)
        p = eru_order.get(L)
        if p:
            target[g] = (gn, ch, L, p)
    LOG(f"[repair] placement references  : {len(target):,} genes at real loci")
    if not target:
        return

    if not (which("diamond") or (which("blastp") and which("makeblastdb"))):
        LOG("[repair] skipped: neither diamond nor blastp/makeblastdb in PATH")
        return

    pep_of = {chrom: pep for _g, _ref, chrom, _c, pep in man}
    tmp = mkdir(os.path.join(wd, "tmp"))
    qfa = os.path.join(tmp, "repair_query.pep")
    sfa = os.path.join(tmp, "repair_target.pep")
    nq = _write_pep(pep_of, query, qfa)
    nt = _write_pep(pep_of, target, sfa)
    LOG(f"[repair] query {nq:,} / target {nt:,} proteins")
    db = os.path.join(tmp, "repairdb")
    bo = os.path.join(tmp, "repair.tsv")
    if not _run_besthit(cfg, sfa, db, qfa, bo, "repair"):
        return
    hits = _read_hits(bo)

    # What each cell already holds.  A move puts a gene into cell (L, chromosome)
    # -- and a cell must hold genes that are all at L *on that chromosome*.
    # Refusing when the cell is taken is the conservative reading of the same
    # principle the rest of the pipeline follows: position authorises, and when
    # position does not separate two candidates, do not merge them.
    cell = defaultdict(list)
    for g, L in gene_loc.items():
        m = gene_meta.get(g)
        if m:
            cell[(L, m[3])].append(g)

    stats = Counter()
    for h, (hgrp, hchrom, hbnd, hiv) in query.items():
        best = None
        for (s, pid, ln, bs) in hits.get(h, ()):
            if pid < MINID or ln < MINLEN:
                continue
            t = target.get(s)
            if not t:
                continue
            tgrp, tchrom, L, (ec, ei) = t
            if tgrp != hgrp:
                stats["refused_different_homology_group"] += 1
                continue
            if tchrom == hchrom:
                stats["refused_same_chromosome_paralog"] += 1
                continue
            lo, hi = hbnd
            if not (lo - WIN <= ei <= hi + WIN):
                stats["refused_outside_the_bracketing_interval"] += 1
                continue
            # the TARGET CELL must be free on h's own chromosome, otherwise the
            # move would put two genes of one chromosome into one cell
            if cell.get((L, hchrom)):
                stats["refused_target_cell_taken_on_this_chromosome"] += 1
                continue
            if best is None or bs > best[1]:
                best = (L, bs, s, pid)
        if best is None:
            continue
        L, _bs, s, pid = best
        old = gene_loc.get(h)
        gene_loc[h] = L
        cell[(L, hchrom)].append(h)
        # works whether h already had a gene_meta entry or was unplaced until now
        gene_meta[h] = ("repaired_to_homoeolog_locus", hiv, hgrp, hchrom,
                        False, hbnd)
        if h in unplaced:
            del unplaced[h]
            stats["repaired_from_unplaced"] += 1
        else:
            stats["repaired_from_synthetic_slot"] += 1
        if old and old != L:
            stats["moved_off_previous_locus"] += 1

    LOG("")
    LOG("[repair] outcome:")
    for k, v in sorted(stats.items()):
        LOG(f"[repair]   {k:48s} {v:>9,}")
    kept = stats["repaired_from_synthetic_slot"] + stats["repaired_from_unplaced"]
    LOG(f"[repair] genes moved onto an allelic locus: {kept:,}")
    cfg.set("_unplaced", unplaced)
    LOG(f"[repair] still unplaced: {len(unplaced):,}")


# --------------------------------------------------------------------------- #
#  STAGE 4 -- refinement
# --------------------------------------------------------------------------- #
def stage_refine(cfg):
    if not cfg.b("refinement", "enabled", True):
        LOG("[refine] disabled in config")
        return
    gene_loc = cfg.get("_gene_loc")
    gene_meta = cfg.get("_gene_meta")
    wd = cfg.wd
    man = load_manifest(wd)

    rows = defaultdict(lambda: defaultdict(list))
    for g, L in gene_loc.items():
        rows[L][gene_meta[g][3]].append(g)
    size = {L: sum(len(v) for v in cd.values()) for L, cd in rows.items()}

    eru_order = cfg.get("_eru_order")
    RICH = cfg.i("refinement", "rich_min", 4)
    SPARSE = cfg.i("refinement", "sparse_max", 2)
    WINDOW = cfg.i("refinement", "window", 10)
    MINID = cfg.f("refinement", "min_identity", 70)
    MINLEN = cfg.i("refinement", "min_alignment", 60)

    def parse_iv(iv):
        if not iv or "|iv|" not in iv:
            return None
        try:
            c, rest = iv.split("|iv|", 1)
            lo, hi = rest.replace("inv", "").split("-")
            return (c, int(lo), int(hi))
        except ValueError:
            return None

    def row_pos(L):
        g0 = next(iter(rows[L].values()))[0]
        src, iv, _gn, _ch, syn, bnd = gene_meta[g0]
        if not syn and L in eru_order:
            c, i = eru_order[L]
            return (c, i, i)
        return parse_iv(iv)

    rich = defaultdict(list)
    for L, sz in size.items():
        if sz >= RICH:
            p = row_pos(L)
            if p:
                rich[p[0]].append((p[1], p[2], L))
    for c in rich:
        rich[c].sort()

    sat, cand_of = [], {}
    for L, sz in size.items():
        if sz > SPARSE:
            continue
        p = row_pos(L)
        if not p:
            continue
        c, lo, hi = p
        cs = [R for (rlo, rhi, R) in rich.get(c, ())
              if rhi >= lo - WINDOW and rlo <= hi + WINDOW]
        if not cs:
            continue
        for _ch, gs in rows[L].items():
            for g in gs:
                sat.append(g)
                cand_of[g] = cs
    LOG(f"[refine] rich rows (>= {RICH}): "
        + ", ".join(f"{c}:{len(v)}" for c, v in sorted(rich.items())))
    LOG(f"[refine] satellite genes considered: {len(sat):,}")
    if not sat:
        return

    if not (which("diamond") or (which("blastp") and which("makeblastdb"))):
        LOG("[refine] skipped: neither diamond nor blastp/makeblastdb in PATH")
        return

    need = set()
    for v in cand_of.values():
        need.update(v)
    db_genes = set()
    for R in need:
        for gs in rows[R].values():
            db_genes.update(gs)
    pep_of = {chrom: pep for _g, _ref, chrom, _c, pep in man}
    tmp = mkdir(os.path.join(wd, "tmp"))
    sfa = os.path.join(tmp, "rich.pep")
    qfa = os.path.join(tmp, "satellite.pep")
    ns = _write_pep(pep_of, db_genes, sfa)
    nq = _write_pep(pep_of, sat, qfa)
    db = os.path.join(tmp, "richdb")
    bo = os.path.join(tmp, "refine.tsv")
    LOG(f"[refine] query {nq:,} / target {ns:,} proteins")
    if not _run_besthit(cfg, sfa, db, qfa, bo, "refine"):
        return
    hits = _read_hits(bo)

    gene2row = {g: L for L, cd in rows.items() for _c, gs in cd.items() for g in gs}
    n = 0
    for g in sat:
        cs = set(cand_of[g])
        best = None
        for s, pid, ln, bs in hits.get(g, ()):
            if pid < MINID or ln < MINLEN:
                continue
            R = gene2row.get(s)
            if R in cs and (best is None or bs > best[1]):
                best = (R, bs)
        if not best or best[0] == gene_loc[g]:
            continue
        oldL = gene_loc[g]
        ch = gene_meta[g][3]
        try:
            rows[oldL][ch].remove(g)
        except (KeyError, ValueError):
            pass
        rows[best[0]][ch].append(g)
        gene_loc[g] = best[0]
        src, iv, gn, _ch, syn, bnd = gene_meta[g]
        gene_meta[g] = ("merged_homology", iv, gn, ch, False, bnd)
        n += 1
    LOG(f"[refine] merged {n:,} satellite genes back into rich rows "
        f"(identity >= {MINID:g}, alignment >= {MINLEN})")


# --------------------------------------------------------------------------- #
#  STAGE 5 -- report / write outputs
# --------------------------------------------------------------------------- #
def write_outputs(cfg):
    wd, od = cfg.wd, cfg.od
    mkdir(od)
    man = load_manifest(wd)
    gene_loc = cfg.get("_gene_loc")
    gene_meta = cfg.get("_gene_meta")
    prefix = cfg.s("output", "prefix", "polyallele")

    xpos, gspan = {}, {}
    for _g, _ref, chrom, _c, _p in man:
        bed = os.path.join(wd, "target", chrom + ".bed")
        for i, (_s, g) in enumerate(read_bed(bed)):
            xpos[g] = i
        with open(bed) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 4:
                    try:
                        gspan[p[3]] = (int(p[1]), int(p[2]))
                    except ValueError:
                        pass

    # Copy identity -- a paralog must be *labelled* as the kind of paralog it is.
    # Distance is measured from the cell's representative, on the same chromosome.
    TANDEM = cfg.f("output", "tandem_kb", 10) * 1000.0
    PROX = cfg.f("output", "proximal_kb", 1000) * 1000.0

    def copy_class(rep, g):
        a, b = gspan.get(rep), gspan.get(g)
        if not a or not b:
            return "unknown"
        d = abs(b[0] - a[0])
        if d <= TANDEM:
            return "tandem_copy"
        if d <= PROX:
            return "proximal_copy"
        return "dispersed_copy"

    eloc = {}
    for bed in sorted(glob.glob(os.path.join(wd, "reference", "*.bed"))):
        with open(bed) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 4:
                    eloc[p[3]] = f"{p[0]}:{p[1]}-{p[2]}"

    rows = defaultdict(lambda: defaultdict(list))
    for g, L in gene_loc.items():
        rows[L][gene_meta[g][3]].append(g)

    # Column index comes from the chromosome's POSITION inside its homology group,
    # not from the last letter of its name.  A name like `Chr01` has no subgenome
    # letter at all, so keying on the name would break on any other naming scheme.
    chrom_col, grp_order = {}, defaultdict(list)
    for gname, _ref, chrom, _c, _p in man:
        grp_order[gname].append(chrom)
    for gname, chroms in grp_order.items():
        for k, c in enumerate(sorted(set(chroms))):
            chrom_col.setdefault(c, k)
    n_cols = max(chrom_col.values(), default=0) + 1
    cols = [f"Allele_{LETTERS[k]}" for k in range(n_cols)]

    def is_synthetic(L):
        return ("|iv|" in L) or ("|same|" in L)

    def col_of(chrom):
        k = chrom_col.get(chrom)
        return cols[k] if k is not None and k < len(cols) else "Extra_Copies"

    split_syn = cfg.b("output", "split_synthetic", True)
    f_long = os.path.join(od, f"{prefix}.allele_long.tsv")
    f_syn = os.path.join(od, f"{prefix}.synthetic_rows.tsv")
    o_long = open(f_long, "w", newline="")
    w_long = csv.writer(o_long, delimiter="\t", lineterminator="\n")
    w_long.writerow(["LocusID", "HomologyGroup", "Ref_Gene", "Ref_Locus",
                     "Subgenome", "Chromosome", "GeneID", "Rank", "CopyClass",
                     "Source", "Interval", "Synthetic"])
    o_syn = open(f_syn, "w", newline="") if split_syn else None
    w_syn = csv.writer(o_syn, delimiter="\t", lineterminator="\n") if o_syn else None
    if w_syn:
        w_syn.writerow(["LocusID", "HomologyGroup", "Interval", "Subgenome",
                        "Chromosome", "GeneID", "Rank"])

    for L in sorted(rows):
        syn = is_synthetic(L)
        rgene = "" if syn else L
        for chrom in sorted(rows[L]):
            gs = sorted(rows[L][chrom],
                        key=lambda x: (0 if gene_meta[x][0] == "anchor" else 1, xpos.get(x, 10**9)))
            gn = gene_meta[gs[0]][2]
            sub = LETTERS[chrom_col.get(chrom, 0)] if chrom in chrom_col else "?"
            for i, g in enumerate(gs):
                src, iv, _gn, _c, sy, _b = gene_meta[g]
                cls = "representative" if i == 0 else copy_class(gs[0], g)
                w_long.writerow([L, gn, rgene if not syn else "NA",
                                 eloc.get(rgene, "NA"), sub, chrom,
                                 g, i + 1, cls, src, iv, int(sy)])
                if w_syn and syn:
                    w_syn.writerow([L, gn, iv, sub, chrom, g, i + 1])
    o_long.close()
    if o_syn:
        o_syn.close()

    # wide matrix
    f_wide = os.path.join(od, f"{prefix}.allele_matrix.tsv")
    n_cells = n_single = n_extra = 0
    class_tally = Counter()
    with open(f_wide, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["LocusID", "HomologyGroup", "Ref_Gene", "Ref_Locus"]
                   + cols + ["Extra_Copies", "Extra_Classes", "NCopies",
                             "Synthetic", "Evidence", "Merged"])
        for L in sorted(rows):
            syn = is_synthetic(L)
            rgene = "" if syn else L
            row = {c: "NA" for c in cols}
            extra, nc, gn = [], 0, None
            classes = Counter()
            allg = [g for chrom in rows[L] for g in rows[L][chrom]]
            n_merged = sum(1 for g in allg
                           if gene_meta[g][0] == "merged_homoeolog")
            # The row's own positional evidence.  A synthetic row has no
            # reference orthologue at all; an anchored row holds at least one
            # collinear anchor; anything else reached a real reference locus by
            # interpolation, which fixes the locus but not the pairing.
            if syn:
                ev = "slot"
            elif any(gene_meta[g][0] == "anchor" for g in allg):
                ev = "anchor"
            else:
                ev = "interval"
            for chrom in sorted(rows[L]):
                gs = sorted(rows[L][chrom],
                            key=lambda x: (0 if gene_meta[x][0] == "anchor" else 1, xpos.get(x, 10**9)))
                gn = gn or gene_meta[gs[0]][2]
                row[col_of(chrom)] = gs[0]
                nc += len(gs)
                n_cells += 1
                if len(gs) == 1:
                    n_single += 1
                else:
                    n_extra += len(gs) - 1
                    for g in gs[1:]:
                        extra.append(f"{chrom}:{g}")
                        c = copy_class(gs[0], g)
                        classes[c] += 1
                        class_tally[c] += 1
            w.writerow([L, gn, rgene if not syn else "NA",
                        eloc.get(rgene, "NA")] + [row[c] for c in cols]
                       + [",".join(extra) if extra else "NA",
                          ";".join(f"{k}={v}" for k, v in sorted(classes.items()))
                          if classes else "NA",
                          nc, int(syn), ev, n_merged])

    # unplaced
    placed = set(gene_loc)
    f_un = os.path.join(od, f"{prefix}.unplaced_genes.tsv")
    n_un = 0
    with open(f_un, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["GeneID", "Chromosome", "GeneOrderInChrom", "Start", "End"])
        for _g, _ref, chrom, _c, _p in man:
            bed = os.path.join(wd, "target", chrom + ".bed")
            with open(bed) as fh:
                for i, line in enumerate(fh):
                    p = line.rstrip("\n").split("\t")
                    if len(p) >= 4 and p[3] not in placed:
                        w.writerow([p[3], chrom, i, p[1], p[2]])
                        n_un += 1

    # stats
    tot = len(xpos)
    anc = sum(1 for g in gene_meta if gene_meta[g][0] == "anchor")
    syn_n = sum(1 for g in gene_meta if gene_meta[g][4])
    h = Counter()
    for L in rows:
        h[sum(len(v) for v in rows[L].values())] += 1
    lines = ["=" * 84,
             f" PolyAllele v{VERSION} -- allele table summary",
             "=" * 84,
             f" target genes                    : {tot:,}",
             f" genes placed on a locus         : {len(placed):,} "
             f"({100.0*len(placed)/max(tot,1):.2f}%)",
             f"   tier-1 anchored               : {anc:,}",
             f"   interpolated, reference locus : {len(placed)-anc-syn_n:,}",
             f"   interpolated, synthetic slot  : {syn_n:,}",
             f" genes NOT placed                : {n_un:,} "
             f"({100.0*n_un/max(tot,1):.2f}%)",
             f" loci (rows)                     : {len(rows):,}",
             f" filled cells                    : {n_cells:,}",
             f" cells with exactly 1 gene       : {n_single:,} "
             f"({100.0*n_single/max(n_cells,1):.2f}%)",
             f" extra copies                    : {n_extra:,}",
             f" inflation (genes / cells)       : {len(placed)/max(n_cells,1):.4f}x",
             "",
             " copy identity (extras, by distance from the cell representative):",
             f"   tandem     (<= {TANDEM/1000:g} kb)  : {class_tally.get('tandem_copy',0):,}",
             f"   proximal   (<= {PROX/1000:g} kb) : {class_tally.get('proximal_copy',0):,}",
             f"   dispersed  (> {PROX/1000:g} kb)   : {class_tally.get('dispersed_copy',0):,}",
             f"   unknown                          : {class_tally.get('unknown',0):,}",
             "",
             " row-size distribution:"]
    for k in sorted(h)[:16]:
        lines.append(f"   {k:4d} gene(s): {h[k]:8,}")
    lines.append(f"   singletons: {h.get(1,0):,} "
                 f"({100.0*h.get(1,0)/max(sum(h.values()),1):.1f}% of rows)")
    txt = "\n".join(lines) + "\n"
    with open(os.path.join(od, f"{prefix}.stats.txt"), "w") as o:
        o.write(txt)
    LOG("")
    LOG(txt)
    LOG(f"[output] {f_wide}")
    LOG(f"[output] {f_long}")
    if split_syn:
        LOG(f"[output] {f_syn}")
    LOG(f"[output] {f_un}")


def stage_merge(cfg):
    """Bring a lone gene home.

    A row holding exactly ONE gene can mean two very different things: the gene
    really is alone in the genome, or its homoeologs were placed in OTHER rows --
    which happens wherever the anchors are dense on one chromosome and sparse on
    the other, so the same ancestral locus is written at two different slots.

    The second case makes a locus look single-copy when it is not, and that is
    exactly what a downstream reader takes "one gene in this row" to mean.  This
    stage moves such a gene into the row that already holds its partner, and is
    allowed to do so only when ALL of these hold:

      * the gene is ALONE in its row;
      * its best non-self BLASTP hit at >= min_identity over >= min_alignment
        amino acids is a gene that IS in the table;
      * same homology group but a DIFFERENT chromosome.  A partner on the gene's
        own chromosome is a paralog, not an allele, and is refused;
      * the destination cell -- (partner's locus, the gene's OWN chromosome) --
        is EMPTY.  Filling an occupied cell would put two genes of one chromosome
        into one cell, which is the over-merge this pipeline exists to avoid, and
        is what made the unbounded v4.2 collision rule fail validation.

    Position authorises; similarity only chooses among what position allows.
    """
    if not cfg.b("merge", "enabled", True):
        LOG("[merge] disabled in config")
        return
    gene_loc = cfg.get("_gene_loc")
    gene_meta = cfg.get("_gene_meta")
    unplaced = cfg.get("_unplaced") or {}
    eru_order = cfg.get("_eru_order") or {}
    if not gene_loc:
        LOG("[merge] no placement state; nothing to do")
        return

    MINID = cfg.f("merge", "min_identity", 90)
    MINLEN = cfg.i("merge", "min_alignment", 80)

    # ---- state ---------------------------------------------------------------
    row_size = Counter(gene_loc.values())
    cell = defaultdict(set)                     # (locus, chromosome) -> {gene}
    for g, L in gene_loc.items():
        m = gene_meta.get(g)
        if m:
            cell[(L, m[3])].add(g)

    # A tier-1 anchor is the STRONGEST positional evidence there is: JCVI put this
    # gene on this reference locus from collinearity itself, not from an
    # interpolation rule.  When such a gene is alone in its row and a >=90%
    # relative sits elsewhere, position and similarity disagree -- and here
    # position wins, because overriding it would mean trusting a BLAST hit over a
    # collinear block.  Anchored genes are therefore never moved.
    #
    # They are still SCREENED, though: the query set is every lone gene, anchored
    # or not, so that changing the anchor policy later does not invalidate the
    # cached alignment and cost another hour of BLAST.
    singles = sorted(g for g, L in gene_loc.items()
                     if row_size[L] == 1 and g not in unplaced)
    movable = [g for g in singles
               if (gene_meta.get(g) or ("",))[0] != "anchor"]
    LOG(f"[merge] genes alone in their row : {len(singles):,} "
        f"of which {len(singles)-len(movable):,} are tier-1 anchors, left where "
        f"collinearity put them")
    if not singles:
        cfg.set("_merged", set())
        return

    # ---- who may receive -----------------------------------------------------
    # every placed gene, including one that is itself alone: two lone genes that
    # are each other's homoeolog belong together, and merging them turns two
    # misleading rows into one honest two-column row.
    target = {}
    for g, L in gene_loc.items():
        m = gene_meta.get(g)
        if not m:
            continue
        src, iv, gn, ch, syn, bnd = m
        target[g] = (gn, ch, L, src)
    LOG(f"[merge] placement references  : {len(target):,} genes at real loci")
    if not target:
        cfg.set("_merged", set())
        return

    if not (which("diamond") or (which("blastp") and which("makeblastdb"))):
        LOG("[merge] skipped: neither diamond nor blastp/makeblastdb in PATH")
        cfg.set("_merged", set())
        return

    man = load_manifest(cfg.wd)
    pep_of = {chrom: pep for _g, _ref, chrom, _c, pep in man}
    tmp = mkdir(os.path.join(cfg.wd, "tmp"))
    qfa = os.path.join(tmp, "merge_query.pep")
    sfa = os.path.join(tmp, "merge_target.pep")
    nq = _write_pep(pep_of, singles, qfa)
    nt = _write_pep(pep_of, target, sfa)
    LOG(f"[merge] query {nq:,} / target {nt:,} proteins")
    db = os.path.join(tmp, "mergedb")
    bo = os.path.join(tmp, "merge.tsv")
    if not _run_besthit(cfg, sfa, db, qfa, bo, "merge"):
        cfg.set("_merged", set())
        return
    hits = _read_hits(bo)

    # ---- the move ------------------------------------------------------------
    stats = Counter()
    for g in singles:
        if (gene_meta.get(g) or ("",))[0] == "anchor":
            stats["left_alone_tier1_anchor"] += 1
            continue
        hgrp = gene_meta.get(g, (None, None, None, None, None, None))[2]
        hchrom = gene_meta.get(g, (None, None, None, None, None, None))[3]
        if not hgrp or not hchrom:
            stats["skipped_no_metadata"] += 1
            continue
        if row_size[gene_loc[g]] != 1:
            stats["skipped_row_already_grew"] += 1
            continue
        best = None
        for (s, pid, ln, bs) in hits.get(g, ()):
            if pid < MINID or ln < MINLEN:
                continue
            t = target.get(s)
            if not t:
                continue
            tgrp, tchrom, L, _tsrc = t
            if tgrp != hgrp:
                stats["refused_different_homology_group"] += 1
                continue
            if tchrom == hchrom:
                stats["refused_same_chromosome_paralog"] += 1
                continue
            if cell.get((L, hchrom)):
                stats["refused_target_cell_taken_on_this_chromosome"] += 1
                continue
            if best is None or bs > best[1]:
                best = (L, bs, s, pid)
        if best is None:
            stats["no_eligible_partner"] += 1
            continue
        L, _bs, _s, _pid = best
        old = gene_loc[g]
        was_synthetic = (gene_meta.get(g) or ("",))[0]
        src, iv, gn, ch, syn, bnd = gene_meta[g]
        cell[(old, hchrom)].discard(g)
        row_size[old] -= 1
        gene_loc[g] = L
        gene_meta[g] = ("merged_homoeolog", iv, gn, ch, False, bnd)
        cell[(L, hchrom)].add(g)
        row_size[L] += 1
        stats["merged_into_an_existing_row"] += 1
        if str(was_synthetic).startswith("interp_synthetic") or \
           str(was_synthetic).endswith("_pushed"):
            stats["merged_from_a_synthetic_slot"] += 1
        if row_size[old] == 0:
            stats["source_row_became_empty"] += 1

    LOG("")
    LOG("[merge] outcome:")
    for k, v in sorted(stats.items()):
        LOG(f"[merge]   {k:48s} {v:>9,}")
    merged = {g for g in singles
              if gene_meta.get(g, ("",))[0] == "merged_homoeolog"}
    cfg.set("_merged", merged)
    LOG(f"[merge] genes moved to the row holding their homoeolog: {len(merged):,}")
    LOG(f"[merge] rows removed : "
        f"{stats['source_row_became_empty']:,}   remaining : "
        f"{sum(1 for _l, n in Counter(gene_loc.values()).items() if n):,}")


def stage_report(cfg):
    write_outputs(cfg)


# --------------------------------------------------------------------------- #
#  main
# --------------------------------------------------------------------------- #
def run_stage(cfg, stage):
    if stage == "prep":
        stage_prep(cfg)
    elif stage == "project":
        stage_project(cfg)
    elif stage == "interp":
        stage_interp(cfg)
    elif stage == "repair":
        stage_repair(cfg)
    elif stage == "refine":
        stage_refine(cfg)
    elif stage == "merge":
        stage_merge(cfg)
    elif stage == "report":
        stage_report(cfg)
    else:
        raise SystemExit(f"unknown stage: {stage}")


def emit_shell(cfg):
    """Print shell assignments for the runtime section, for run_polyallele.sh."""
    def q(s):
        return "'" + str(s).replace("'", "'\\''") + "'"
    kv = {
        "POLYALLELE_CONDA_ENV": cfg.s("runtime", "conda_env", ""),
        "POLYALLELE_SCHEDULER": cfg.s("runtime", "scheduler", "none").lower(),
        "POLYALLELE_QUEUE": cfg.s("runtime", "pbs_queue", "comput"),
        "POLYALLELE_PPN": str(cfg.i("runtime", "pbs_ppn", 8)),
        "POLYALLELE_WALLTIME": cfg.s("runtime", "pbs_walltime", "12:00:00"),
        "POLYALLELE_JOBNAME": cfg.s("runtime", "job_name", "polyallele"),
        "POLYALLELE_WORKDIR": cfg.p("output", "work_dir", "./polyallele_work"),
        # the tables that sit NEXT TO the matrix are built by bin/extra_tables.py
        "POLYALLELE_OD": cfg.p("output", "out_dir", "./polyallele_out"),
        "POLYALLELE_PREFIX": cfg.s("output", "prefix", "polyallele"),
        "POLYALLELE_EXTRAS": "1" if cfg.b("copyclass", "enabled", True) else "0",
        "POLYALLELE_ALLVALL": cfg.p("copyclass", "allvall", "") or "",
        "POLYALLELE_IDENTITY": str(cfg.f("copyclass", "min_identity", 90)),
        "POLYALLELE_COMPID": str(cfg.f("copyclass", "competition_identity", 95)),
        "POLYALLELE_MINLEN": str(cfg.i("copyclass", "min_alignment", 60)),
    }
    for k, v in kv.items():
        print(f"{k}={q(v)}")


def main():
    ap = argparse.ArgumentParser(
        description=f"PolyAllele v{VERSION} -- position-first allele matrix builder")
    ap.add_argument("--config", default="config.ini", help="path to config.ini")
    ap.add_argument("--emit-shell", action="store_true",
                    help="print runtime.* as shell assignments and exit "
                         "(used by run_polyallele.sh)")
    ap.add_argument("stage", nargs="?", default="all",
                    choices=["prep", "project", "interp", "repair", "refine",
                             "merge", "report", "all"],
                    help="stage to run; prerequisites are run automatically "
                         "(default: all)")
    args = ap.parse_args()

    cfg = Cfg(args.config)
    if args.emit_shell:
        emit_shell(cfg)
        return

    cfg.wd = mkdir(cfg.p("output", "work_dir", "./polyallele_work"))
    cfg.od = mkdir(cfg.p("output", "out_dir", "./polyallele_out"))

    global LOG
    LOG = Log(os.path.join(cfg.od, "polyallele.log"))
    LOG("")
    LOG("=" * 84)
    LOG(f" PolyAllele v{VERSION}   config={cfg.path}")
    LOG(f" work_dir = {cfg.wd}")
    LOG(f" out_dir  = {cfg.od}")
    LOG("=" * 84)

    PIPE = ["prep", "project", "interp", "repair", "refine", "merge", "report"]
    run = PIPE if args.stage == "all" else PIPE[:PIPE.index(args.stage) + 1]

    # Rebuild the state a resumed run needs.  `prep` and `project` are both
    # idempotent, so re-running the prefix is cheap and always safe.
    refchroms = sorted(os.path.basename(b)[:-4]
                       for b in glob.glob(os.path.join(cfg.wd, "reference", "*.bed")))
    have_prep = bool(refchroms) and os.path.exists(os.path.join(cfg.wd, "manifest.tsv"))
    if have_prep:
        cfg.set("_refchroms", refchroms)
        cfg.set("_targetchroms", sorted(os.path.basename(f)[:-4]
                                        for f in glob.glob(os.path.join(cfg.wd, "target", "*.cds"))))
        cfg.set("_groups", resolve_groups(cfg))
        eru_order, eru_lists = load_reference_order(cfg.wd)
        cfg.set("_eru_order", eru_order)
        cfg.set("_eru_lists", eru_lists)

    for s in run:
        if s == "prep" and have_prep and args.stage == "all":
            LOG("")
            LOG("[skip] prep -- reference frames and manifest already present")
            continue
        LOG("")
        LOG(f"--- stage {s} ---")
        t0 = time.time()
        run_stage(cfg, s)
        LOG(f"--- stage {s} done in {int(time.time()-t0)}s ---")

    LOG("")
    LOG("[done]")
    LOG.close()


if __name__ == "__main__":
    main()
