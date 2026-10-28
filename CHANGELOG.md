# Changelog

## v4.4

### Changed — two groups are written, because one cannot serve both purposes

A gene can be **alone at its locus and still have a `>=90%` relative at another
locus**. Such a gene is not single copy, and its row is nonetheless correct. That
is uncomfortable, because a quantification tool reads the rows of a cluster file
and treats a one-gene row as single copy: two genes that are really brothers but
sit in separate one-gene rows do **not** compete, and their expression comes out
inflated.

Merging loci until brothers share a row fixes the inflation and destroys the locus
structure — in a genome where 46% of genes have a 100%-identical relative,
chaining took the largest group to **30,656 genes**. Measured on the reference
dataset, the two groupings disagree for **97% of one-gene rows**:

| | no relative | has a relative |
|---|---|---|
| locus holds 1 gene | 533 | **17,807** |
| locus holds 2+ | 596 | 317,275 |

So both are now written, and neither is asked to do the other's job:

| file | cluster is | for |
|---|---|---|
| `<prefix>.locus2gene.tsv` | an ancestral locus (position) | comparative work |
| `<prefix>.competition2gene.tsv` | genes whose reads cannot be told apart | **quantification** |
| `<prefix>.gene2copyclass.tsv` | both, plus `RowIsSingleCopy` | reading |

`RowIsSingleCopy` is the flag that says what the old table could not: **1 only
when the locus holds one gene AND that gene has no relative at all.**

### Added — the competition threshold is measured, not assumed

`copyclass.competition_identity` is the identity above which two genes fight over
one read — a different question from "are these duplicates", and higher, because
at 150 bp a 90% pair leaves about fifteen mismatches and an aligner resolves that:

| T | components | genes alone | largest cluster |
|---|---|---|---|
| 90% | 51,327 | 16,454 | 20,414 |
| **95%** | 75,447 | 33,554 | 4,486 |
| 98% | 116,689 | 66,452 | 204 |
| 99% | 158,725 | 102,337 | 100 |

Lower means more genes compete — less inflation, more over-grouping. Pick it to
match what the quantifier can tell apart.

### Fixed — `CopyClass` and the locus size are independent, and now say so

The v4.3 extra tables derived `Role` from `Source` and had no way to express the
difference between "one gene in this row" and "a single-copy gene". Both are fixed,
and `test/test_extra_tables.py` pins down each case, including that a gene alone at
its locus but with a relative elsewhere is **not** flagged.

### Fixed — every output file was written with CRLF line endings

`csv.writer` defaults to `lineterminator="\r\n"`, so every TSV the pipeline and
the extra tables produced had Windows line endings. On Linux that breaks
`awk '{print $NF}'`, `cut`, and any `line.split("\t")[-1]`: the last field keeps a
stray `\r`.

It broke this project's own profiler, which read a column header as `'Source\r'`,
failed to find `Source`, and reported every value in that column as missing —
silently, as a `?`. That is the failure mode to watch for: not a crash, a
plausible-looking wrong answer.

Every `csv.writer` now sets `lineterminator="\n"`, and `test/test_formats.py`
asserts both that each of them says so and that the result really is LF. A test
that only checked one file would not have caught the eight call sites.

### Fixed — the copy-class table must be built from BLASTP, not DIAMOND

The first build used DIAMOND `--sensitive` for the genome-wide screen. Checking
its own answer against BLASTP showed the answer was wrong:

```
  DIAMOND said single_copy        : 22,609
  BLASTP finds a >=90% relative   :  5,558  (24.58%)   <- FALSE single-copy calls
  BLASTP agrees: truly alone      : 17,051
```

`--sensitive` is not sensitive enough for a question whose entire content is
"nothing else in this genome is like this gene". The recipe in `config.ini` now
specifies BLASTP, and records why.

The same check exposed a second, larger point about the data:

```
  alone within their homology group       : 42,659
  of those, with a cross-group relative   : 24,045  (56%)
```

**56% of the genes that are alone inside their homology group do have a `>=90%`
relative in another group.** Two consequences: screening only within groups would
have called all 42,659 single copy, and "different homology group" cannot be
treated as automatically impossible when adjudicating a merge — it is strong
evidence, not proof. The definitive table ends with **18,614 single_copy** and
**317,597 multi_copy** genes.

### Fixed — two defects in the extra tables

* the comment line of `unplaced_matrix.tsv` was written as two comma-separated
  fields, so it came out as two TAB-separated columns and read back as part of
  the header;
* `Role` in `locus2gene.tsv` was derived from `Source` instead of `Rank`, which
  left it blank for every interpolated gene — that is, most of them.

### Added — `merge`: a lone gene joins the row that holds its homoeolog

A row with one gene gets read as "nothing else in the genome is this gene's
allele", and it is easy to believe that means the gene is single copy. Usually it
does not. Measuring the table showed **4,638 of the 22,690 one-gene rows** exist
only because the two chromosomes' anchors are dense on one side and sparse on the
other, so one ancestral locus was written at two different slots.

The new `merge` stage moves such a gene into the row that already holds its
partner, and refuses unless **all** of these hold:

| gate | why |
|---|---|
| the gene is alone in its row | otherwise it is not a lone gene |
| its best non-self hit (>= `min_identity`, >= `min_alignment` aa) is in the table | there must be somewhere to go |
| same homology group, **different chromosome** | a partner on the same chromosome is a paralog, not an allele |
| destination cell `(partner's locus, gene's own chromosome)` is **empty** | filling it would put two genes of one chromosome in a cell -- the over-merge that made v4.2's unbounded collision rule fail validation |

`merge.enabled = false` reproduces the previous behaviour exactly.

**Measured on the reference dataset** (sugarcane ROC22 x Erianthus rufipilus,
336,211 genes):

| | before | after |
|---|---|---|
| rows | 57,547 | **54,077** |
| one-gene rows | 22,690 | **17,674** |
| one-gene rows whose gene has no `>=90%` relative | 52.4% | **67.3%** |
| one-gene rows whose gene's allele was placed in another row | 17.4% | **8.4%** |
| genes placed / cells / inflation | unchanged: 312,639 / 261,372 / 1.1961x | |

Gates fired as intended: 114,920 refusals for a different homology group, 24,645
for a same-chromosome partner, 7,640 because the destination cell was taken, and
330 lone **tier-1 anchors** left where collinearity put them — a collinear block
outranks a BLAST hit, so an anchored gene is never moved.

External validation, unchanged or better on every axis:

| check | v4.1 | v4.3 |
|---|---|---|
| curated WOX family, over-merge | 0 | **0** |
| cross-group merges | 0 | **0** |
| agreement with an independent method, our grouping reproduced | 81.44% | 78.33% |
| our disputed merges that are provably impossible | 6.54% | 9.76% |
| 1:1 rate / inflation | ~84% / ~1.19x | **85.03% / 1.1961x** |

The agreement figure is the one number that is lower than v4.1, by 3.1 points on
a 300,000-pair sample. It is reported rather than hidden: the reference method is
itself 34.52% provably-impossible by the same criterion, mostly by merging across
homology groups, which this pipeline never does.

### Added — `Evidence` column, and three tables beside the matrix

The matrix states which ancestral locus a gene belongs to. It cannot state whether
the gene is single copy in the genome, because those are different questions. New:

| file | contents |
|---|---|
| `Evidence` column in the matrix | `anchor` / `interval` / `slot` — how strong the positional evidence is, plus `Merged` |
| `<prefix>.unplaced_matrix.tsv` | genes that entered no row, arranged as a matrix too: rows are >=90% components over the unplaced genes only, columns are chromosomes |
| `<prefix>.locus2gene.tsv` | `LocusID <TAB> GeneID`, table genes only |
| `<prefix>.gene2copyclass.tsv` | **every** gene, `single_copy` or `multi_copy`, from similarity alone |

The last one deliberately ignores position, so its clusters are **similarity**
components, not allele groups — one cluster can span several ancestral loci, and
the file header says so. Built by `bin/extra_tables.py`, which only reads what the
pipeline wrote, so a problem there cannot damage the matrix.

### Added — input format compatibility

Every reader now normalises gene IDs through `input.gene_id_strip`, so one gene
written two ways still matches. A gene ID has to agree between the target BED and
the JCVI anchors or **nothing** is anchored, and the failure then surfaces much
later as an empty result rather than as an error.

| reader | now accepts |
|---|---|
| GFF/GTF | `gene`, `mRNA`, `transcript`, `CDS`; `ID=`, `gene_id=`, `Name=`, `gene=`, `locus_tag=`, `Parent=`; quoted values; CRLF |
| BED | tab or space separated; `track` / `browser` / `#` lines; extra columns; CRLF |
| FASTA | legacy `;` comments; blank lines; CRLF; description text after the ID |
| anchors | `#` comments; missing score column |
| config | a UTF-8 BOM (a config written on Windows) |

`prep` now fails immediately, with an actionable message, when no target
chromosomes are found — instead of producing an empty matrix forty minutes later.

### Fixed — `run_polyallele.sh` under `set -u`

`conda activate` runs `activate.d` scripts that reference variables which may be
unset, and under `set -u` that aborts the job. The activation is wrapped; a helper
script that located its aligner through conda died in twelve seconds with an empty
log because of exactly this.

## v4.2

### Changed — `interp_distributed` is off by default (it caused half the over-merging)

Two code changes failed to move the provably-bad merge rate (15.12% → 14.14% →
14.06%), so instead of a third guess the disagreements were **attributed back to
the stage that produced them**. The long table records a `Source` per gene, so
the question "which stage makes the bad merges?" can be answered by measurement
rather than by another two-hour run:

```
  interp_distributed_copy + interp_distributed_copy      3,964   26.6%
  interp_distributed      + interp_distributed_copy      2,005   13.5%
  interp_distributed_inv_copy + interp_distributed_inv_copy  450   3.0%
  copy_at_anchored_locus  + copy_at_anchored_locus       1,606   10.8%
  anchor                  + merged_homology              1,493   10.0%
  merged_homology         + merged_homology                609    4.1%
  anchor                  + anchor                         326    2.2%
```

Counted by source, `interp_distributed_copy` appears in 10,663 bad merges,
more than any other. That branch was ADDED in v4.2: when an interval holds more
target genes than the reference offers loci, it spread them over the loci and
called the surplus extra copies. The reasoning was that this avoids manufacturing
singleton rows. The reasoning was wrong:

> **When the reference offers four loci and the target has ten genes, calling six
> of them "copies at those loci" is an inference, not evidence -- and the position
> of those six is exactly what is unknown.**

`interpolation.distribute_surplus = false` is now the default: the first `g_gap`
genes take the loci in order and the rest get their own slot. The old behaviour
is still available by setting it to `true`.

**The method is the point.** Three rounds were spent fixing causes that were
plausible rather than measured. Attribution over an existing column answered the
question in two minutes. When several stages could explain a symptom, measure the
contributions before changing any of them.

### Fixed — `repair` had no gate on the TARGET CELL

The three original gates all constrain the **partner**: same homology group,
different chromosome, positionally compatible. None of them looks at where the
gene is going.

That is a hole. A cell is `(locus, chromosome)`, and moving a gene into it is
only sound if the cell is free **on the gene's own chromosome**. Otherwise the
move puts two genes of one chromosome into one cell -- over-merging, which is
the fatal kind.

```
[repair]   refused_target_cell_taken_on_this_chromosome         ...
```

This is what the third-party comparison was pointing at. Bounding the collision
rule by gene order (below) only moved the provably-bad merge rate from 15.12% to
14.14% -- barely at all -- while `repair` was moving *more* genes, not fewer.

**The general lesson, which took three attempts to learn here:** when a stage has
several gates, check that between them they cover **every** way the operation can
be wrong -- not just the ways that are easy to state.

### Changed — the collision rule is bounded by gene order

v4.2 originally absorbed **every** same-chromosome collision into the cell. That
was wrong, and the third-party comparison caught it: intervals in reference-locus
space can overlap when the projection is non-monotonic, so two genes far apart
along a chromosome can "collide" without being at the same locus.

|                              | v4.1 | v4.2 (unbounded) | v4.2 (bounded) |
|------------------------------|------|------------------|----------------|
| provably-bad merges          | 6.54%| **15.12%**       | measured below |
| third-party agreement (v4)   |81.44%| **66.47%**       | measured below |
| third-party agreement (them) |75.22%| 76.13%           |                |

A collision is now absorbed only when the gene is within
`interpolation.tandem_gene_window` (default 5) gene models of one already in that
cell -- which is what a real tandem array looks like. Everything else gets its
own slot, as before v4.2.

**This is the entry to read if you are wondering whether the curated-family
benchmark is enough.** It is not. On the WOX family v4.2 was a clean success:
over-merge 0, over-split 1.25, identical to v4.1. The WOX loci sit on different
chromosomes, so they cannot exercise the same-chromosome collision rule at all.
A benchmark only tests what its own structure can express.

### Fixed — the screening cache ignored the aligner parameters

`_cache_valid` compared only the query and target **file sizes**, so a run
configured for `max_target_seqs = 20` silently reused a cached `k = 50` result
from an earlier run -- a different set of hits, with no warning. The cache key is
now a digest of the sizes **plus** `evalue`, `max_target_seqs` (per stage) and the
aligner. `test/test_sharding.py` asserts that changing either the query or
`max_target_seqs` changes the key.

### Fixed — `KeyError` in `stage_repair` when repairing an unplaced gene

A gene that entered no row at all has **no `gene_meta` entry** -- it never got
one. `stage_repair` nevertheless offered those genes for repair (they are the
"interval too wide" ones, and repairing them is the point), then wrote
`gene_meta[h][1]` for a key that did not exist:

```
File "polyallele.py", line 1208, in stage_repair
    gene_meta[h] = ("repaired_to_homoeolog_locus", gene_meta[h][1], hgrp, ...
KeyError: 'ROC_gene232936'
```

The candidate tuple now carries the interval, so a gene from either source
produces a complete entry. `test/test_repair_unplaced.py` drives exactly this
path with a mocked aligner and asserts that an unplaced gene is moved onto its
allelic locus, gains a `gene_meta` entry, and leaves the unplaced set.

### Fixed — an unrequested change to `max_target_seqs`

Refactoring the three screens into one helper silently changed the value from
**20 to 50** for `interp` and `refine`. It showed up as the similarity screen
writing 44.9 MB instead of 22.7 MB and reassigning 7,403 genes instead of 8,148 --
different numbers from a supposedly unchanged stage. 20 is the value the whole
pipeline was validated at, so it is restored and the setting is now resolved
**per stage**:

```ini
max_target_seqs        = 20     ; interp and refine -- the validated value
max_target_seqs_repair = 50     ; repair searches the whole placed set
```

**The lesson kept in this entry**: a refactor that changes a default is not a
refactor. The only reason it was caught is that the stage printed a byte count
and a reassignment count that could be compared with the previous run.

### Fixed — the query splitter put every record in the first shard

`_split_fasta` chose a shard from the byte offset at which a record starts, but
accumulated only the length of the **header** lines while dividing by the size of
the **whole file**. On a real proteome (short headers, long sequences) the ratio
crawls, so shard 0 received everything and shards 1..n-1 were created empty. The
screening still produced correct results -- the shards are concatenated -- but
n-1 aligner processes sat idle while one ran at a fraction of the requested
threads, making the stage ~2.4x *slower* than not sharding at all.

It was caught in the first real run, from the file listing:

```
.interp.shard0.pep   12,781,820 bytes
.interp.shard1.pep            0
...
.interp.shard7.pep            0
```

**Why the unit test missed it** is the part worth keeping: the fake query used
4-character sequences, so headers dominated the byte count and a header-only
accumulator looked correct. **A test whose data does not resemble the real data
tests nothing.** `test/test_sharding.py` now writes protein-like records (short
header, 12 lines of sequence) and asserts that the shards come out even, at
n = 2, 4 and 8:

```
n=2  sizes=[32039, 27661]                                     even=True
n=4  sizes=[16119, 15920, 15920, 11741]                       even=True
n=8  sizes=[8159, 7960, 7960, 7960, 7960, 7960, 7960, 3781]   even=True
```

### Fixed — three defects in the `report` stage, found the first time it ever ran

The `report` stage had never been executed end to end before this release, and it
did not survive the first attempt.

- **`KeyError: 'Allele_A'`** — the wide-matrix row was initialised as
  `{f"Allele_{c}": "NA" for c in cols}` where `cols` already held `Allele_A`, so
  the dict keys came out as `Allele_Allele_A`. Any run reaching `report` crashed.
- **Synthetic rows were mis-detected.** `syn = "|iv|" in L` missed the synthetic
  loci built from the same-locus rule (`Chr01|same|7|c0`), which would then be
  written as if they carried a reference gene.
- **The subgenome column came from the last character of the chromosome name**
  (`chrom[-1]`). That works for `Chr1A` and breaks for `Chr01`. The column is now
  taken from the chromosome's **position inside its homology group**, so any
  naming scheme works.

### Added

- **Screening results are cached.** The best-hit screens are by far the most
  expensive part of a run, so their output is kept next to a sidecar recording
  the query and target sizes. A re-run after a crash, a walltime kill, or a
  change to a later stage reuses it; any change to the inputs invalidates it.
  `alignment.force = true` redoes the screen anyway.
- **`alignment` section** (`evalue`, `max_target_seqs`) and **`runtime.aligner`**
  (`blastp` default, `diamond` optional).
- Aligner failures now print their stderr tail instead of failing silently.

### Measured — DIAMOND is **not** a drop-in for blastp here

It is two orders of magnitude faster and was very tempting. On 2,000 queries
against the same candidate set:

| | value |
|---|---|
| best partner identical to blastp | **58.65%** |
| pairs passing id>=70 / len>=60, Jaccard | 86.16% |
| pairs passing id>=90 / len>=80, Jaccard | 85.72% |

In a polyploid many homoeologs are near-identical, so top hits are near-ties and
the two aligners break them differently. The default is therefore `blastp`;
`diamond` is offered for exploration only, and the numbers above are documented
in `config.ini` so nobody has to rediscover them.

### Two fixes to over-splitting, both found by measuring what "single copy" means

### The measurement that motivated this release

A gene the table calls single-copy means *no collinear partner was found*. It does
**not** mean its reads map uniquely there. Measuring that gap against the
intra-group all-vs-all BLASTP:

| what the table calls it | genes | has a >=90% relative in the same group |
|---|---|---|
| single-copy row | 25,027 | **49.0%** |
| lone cell in a multi-copy row | 240,108 | 92.8% |
| multi-copy cell | 45,410 | 95.6% |

(The 92.8% / 95.6% figures are mostly *homoeology* -- those genes sit in rows that
already contain their alleles. The 49.0% is the meaningful number.) Splitting that
49% by where the relative sits:

| the single-copy gene's >=90% relative is... | genes | share |
|---|---|---|
| nowhere (genuinely unique) | 12,761 | 51.0% |
| on the **same** chromosome -- a true paralog | 2,251 | 9.0% |
| on **another chromosome of the same group** -- an **allele** | 10,015 | 40.0% |
| in another group | 0 | 0.0% |

and for those 10,015 the partner was **always** in a synthetic row (77%) or
unplaced (23%) -- never at a placed reference locus. That is over-splitting
produced by the projection, not a biological finding.

Downstream this matters a lot: LineageQuant turns a one-gene row into a size-1
cluster and treats it as an ORPHAN, structurally different from a multi-copy
cluster that shares reads through the EM. So a spuriously split gene stops
competing for its own reads.

### Fixed

- **Same-chromosome collisions are kept in the cell instead of exiled.**
  When two genes on the same chromosome both project to reference locus L, that
  is a genuine tandem/segmental copy at L. The surplus used to be pushed to its
  own synthetic row, which turned a multi-copy locus into a one-gene row. It is
  now recorded as a further copy in the same cell (`Extra_Copies`), capped by
  `interpolation.max_copies_at_anchored_locus`. Set
  `interpolation.keep_same_chromosome_collisions = false` for the old behaviour.

- **New `repair` stage: similarity-guided locus repair.**
  A gene in a synthetic slot, or in an interval too wide to interpolate, is moved
  onto an allelic locus when **all** of these hold:

  * the partner sits at a **real** reference locus (not a synthetic slot);
  * **same homology group, different chromosome** -- homoeologous chromosomes are
    allelic by default. A same-chromosome partner would be a **paralog**, and
    that case is refused outright;
  * BLASTP identity >= `repair.min_identity` (90) over >= `repair.min_alignment`
    (80) amino acids;
  * the partner's locus is **positionally compatible**: it lies inside the
    interval the gene was bracketed by, plus `repair.window` (1) loci of slack.

  Position authorises; similarity only chooses among what position allows. The
  refusals are counted and logged, so the guard rails are auditable:

  ```
  [repair]   refused_same_chromosome_paralog          ...
  [repair]   refused_outside_the_bracketing_interval  ...
  [repair]   refused_different_homology_group         ...
  ```

### Also

- `interp` now reports the genes it could not place, with the interval each was
  bracketed by, so `repair` can act on them.

## v4.1

### Fixed

- **Copies at an anchored locus were being dropped.** When both flanking anchors
  sat on the *same* reference locus, the genes between them were skipped
  (`continue`) instead of being placed. That is biologically the *clearest* case
  there is: two anchors agreeing on locus L, with genes physically between them,
  is direct evidence that those genes are additional copies **at L**. They are now
  assigned to L, appear in `Extra_Copies` with `Rank > 1`, and carry
  `Source = copy_at_anchored_locus`. A cap
  (`interpolation.max_copies_at_anchored_locus`, default 50) stops a runaway
  array from inflating a single cell; any surplus gets its own interval slot.

- **`set -u` conflicted with conda activation.** `run_polyallele.sh` used
  `set -euo pipefail` and then sourced conda's `activate.d` scripts, which
  reference unset variables (`ADDR2LINE: unbound variable`) under nounset. Nounset
  is now relaxed around the activation only. This affected every user on a
  conda-managed cluster.

### Added

- **Copy identity labelling.** A paralog is now reported as the *kind* of paralog
  it is, by physical distance from the cell's representative:

  | `CopyClass` | distance | interpretation |
  |---|---|---|
  | `tandem_copy` | ≤ `output.tandem_kb` (10 kb) | tandem array |
  | `proximal_copy` | ≤ `output.proximal_kb` (1000 kb) | local / segmental duplication |
  | `dispersed_copy` | > `output.proximal_kb` | dispersed duplication |

  Emitted as `CopyClass` per gene in the long table and summarised as
  `Extra_Classes` per row in the wide matrix; the tally is also printed in
  `*.stats.txt`.

- **`input.anchors_dir`** — reuse already-computed JCVI `.anchors` files
  (searched recursively) instead of recomputing the pairwise projection.

- **`interpolation.max_copies_at_anchored_locus`** — cap on how many genes one
  anchored locus may absorb.

- `prep` now derives missing BED files from `input.target_gff`, de-duplicates
  chromosomes when `target_prep_dir` contains more than one `.cds` per
  chromosome (e.g. symlinked JCVI working directories), and looks for `.bed` /
  `.pep` in sibling directories rather than only next to the `.cds`.

## v4.0

First release of the position-first rewrite. Compared with the earlier
similarity-clustering implementations (v1.x) the grouping criterion changes,
which is a breaking change to the output contract.

### Grouping

- **Rows are now ancestral loci of a diploid reference genome**, established by
  1:1 collinear projection. Previously rows were connected components of a
  similarity/synteny graph, closed transitively; those components chain, and a
  single spurious anchor could merge unrelated loci.
- **Columns are subgenomes.** A `(row, column)` cell holds at most one
  representative gene; additional copies at the same locus on the same chromosome
  are recorded explicitly in `Extra_Copies` / `Rank`.
- **The 1:1 property is guaranteed by construction.** A candidate locus already
  occupied by a chromosome is never overwritten.

### New stages

- **`interp` — order-preserving interpolation.** Genes between two collinear
  anchors are placed onto the reference loci inside the interval by order alone.
  Handles **inversions** (reverse walk) and creates **synthetic interval slots**
  when the reference has no room, so the 1:1 property survives regions the
  reference does not cover.
- **`refine` — position-constrained homology refinement.** A gene alone in a small
  row is merged into a nearby rich row only when the row is within `window`
  reference loci **on the same reference chromosome** and BLASTP supports the
  merge. Sequence similarity is never allowed to widen the window; it only decides
  inside it. This repairs over-splitting without reintroducing over-merging.
- **`similarity_first`** optionally resolves `interp_proportional` placements
  inside their interval with BLASTP.

### Output contract

- Wide matrix: `LocusID, HomologyGroup, Ref_Gene, Ref_Locus, Allele_*, Extra_Copies,
  NCopies, Synthetic`.
- Long table gains `Source` (how the locus was assigned) and `Interval`.
- **Synthetic rows are separated by default** (`output.split_synthetic`): they are
  positional annotations with no reference orthologue and should not be presented
  as verified loci.
- New `bin/qc_matrix.py` reports the metrics that decide usability: 1:1 rate,
  inflation, subgenome occupancy, over-merge check, gene-uniqueness check and span
  purity of multi-gene cells.

### Configuration and packaging

- Single `config.ini` drives inputs, chromosome grouping, algorithm parameters,
  outputs, and how to launch (direct / PBS / SLURM).
- One command: `bash run_polyallele.sh`.
- Every stage is idempotent and a later stage runs its prerequisites, so the
  pipeline resumes cleanly after a walltime kill.
- No hardcoded paths, accounts or reference organisms.

### Validation

Validated on sugarcane ROC22 (97 chromosomes in 10 homology groups) against the
diploid *Erianthus rufipilus*:

| metric | this pipeline | similarity-clustering baseline |
|---|---|---|
| multi-gene cells | 11% | 31% |
| inflation (genes/cell) | 1.17× | 1.71× |
| genes in a row of >20 | 0 | 510 |
| largest row | 26 genes | 497 genes |
| genes placed twice | 0 | 0 |
| over-merge, curated WOX family | **0** | 2 |
| within-group over-split, curated WOX family | 1.42 | 1.25 |

Gene placement was additionally checked by hold-out: hiding 20% of the anchors and
re-interpolating recovered the correct locus for 80% of forced 1:1 placements and
put the true locus inside the stated interval for 85% of synthetic slots.

## v1.3 and earlier

Similarity/synteny-graph implementations. Superseded; retained for reference only.
Their output is a gene-family table rather than an allele table, and rows can
contain genes that belong to different ancestral loci.
