# PolyAllele v4

**A position-first pipeline for building standardised allele matrices in polyploid genomes.**

Input: a polyploid genome annotation plus a well-assembled **diploid** relative.
Output: one table where **each row is one ancestral locus** and **each column is one
subgenome copy** — plus an explicit, auditable record of every extra copy and every
gene that could not be placed.

One command:

```bash
bash run_polyallele.sh
```

Everything else lives in `config.ini`.

---

## Table of contents

- [Why this exists](#why-this-exists)
- [What it produces](#what-it-produces)
- [Installation](#installation)
- [Quick start](#quick-start)
- [The configuration file](#the-configuration-file)
- [How the algorithm works](#how-the-algorithm-works)
- [Output files](#output-files)
- [How to read the output](#how-to-read-the-output)
- [Quality control](#quality-control)
- [Parameters that matter](#parameters-that-matter)
- [Limitations — read this before publishing](#limitations--read-this-before-publishing)
- [Resuming after an interruption](#resuming-after-an-interruption)
- [Citation](#citation)
- [License](#license)

---

## Why this exists

A polyploid genome contains several homoeologous copies of most genes, spread over
several chromosomes. To quantify them — and above all to tell *subgenome* expression
apart from *paralog* expression — you need to know which copies belong together.

There are two ways to decide that, and they are **not** equivalent:

| criterion | what it asks | failure mode |
|---|---|---|
| **sequence similarity** | "are these two genes alike?" | recent paralogs are more similar to each other than homoeologs are to each other → **over-merging**: one row silently mixes different genes, and any expression number computed from it is uninterpretable |
| **position (synteny)** | "do these two genes sit at the same place on homoeologous chromosomes?" | a gene whose projection is slightly off lands in the wrong row → **over-splitting**: one gene's copies end up in two rows |

Over-merging is **fatal** — the resulting number has no biological meaning.
Over-splitting is **recoverable** — summing two rows downstream costs resolution,
not correctness.

PolyAllele is therefore built the other way round from similarity-clustering
pipelines:

> **Row identity is decided by collinear position against a diploid reference.
> Sequence similarity is used only as a tie-break inside a window that position has
> already authorised.**

This gives a matrix whose rows are *defined*, whose columns are *subgenomes*, and
whose errors are the recoverable kind.

---

## What it produces

```
LocusID              HomologyGroup  Ref_Gene            Ref_Locus
Erufi.01G0000030     Group_01       Erufi.01G0000030    Chr01:122452-126603
    Allele_A            Allele_B            ...   Allele_L
    ROC_gene107904      ROC_gene075683      ...   NA
    Extra_Copies                  Extra_Classes                    NCopies  Synthetic
    Chr1I:ROC_gene324652          proximal_copy=1;dispersed_copy=2  13       0
```

- **one row = one ancestor locus of the reference genome**
- **one column = one subgenome copy** (the `Allele_A … Allele_L` columns)
- **`Extra_Copies`** = further copies at the *same* locus on the *same* chromosome
  (real duplications; collapse them into the cell downstream)
- **`Extra_Classes`** = what *kind* of duplication each extra is, by distance from
  the cell's representative — see [copy identity](#copy-identity-what-kind-of-paralog)
- **`NCopies`** = how many target genes the row holds in total
- **`Synthetic` = 1** marks rows with **no reference orthologue**: the gene sits in a
  reference-bracketed interval, so its position is known but its ancestor locus is not

The pipeline also emits a long table (one line per gene, with `Rank`, `CopyClass`
and `Source`) and a list of genes it could not place at all.

---

## Installation

```bash
conda env create -f environment.yml
conda activate polyallele
```

The environment provides `jcvi`, `last`, `blast` and `gffread`.
If you do not use conda, install those four tools yourself; the pipeline calls them
by name from `PATH`.

Check the installation:

```bash
python bin/polyallele.py --help
bash run_polyallele.sh config.ini prep      # runs stage 1 only
```

See [INSTALL.md](INSTALL.md) for a manual route and for cluster notes.

---

## Quick start

1. **Edit `config.ini`.** At minimum set the reference genome, the target genome,
   and how chromosomes are grouped (see below).
2. **Run:**

```bash
bash run_polyallele.sh
```

3. **Inspect:**

```bash
python bin/qc_matrix.py --matrix polyallele_out/polyallele.allele_matrix.tsv \
                        --bed-dir  polyallele_work/target \
                        --report   polyallele_out/qc.txt
```

On an HPC, set `[runtime] scheduler = pbs` (or `slurm`) and the same command will
write a job script and submit it for you.

---

## The configuration file

`config.ini` has eight sections. The ones you must touch:

### `[input]`
Either give the whole-genome files and let the pipeline extract CDS with `gffread`:

```ini
target_gff   = /data/XTT.gff3
target_fasta = /data/XTT.fasta
```

or point at a directory of **pre-extracted per-chromosome files**
(`<chrom>.cds`, `<chrom>.bed`, optionally `<chrom>.pep`) — faster and reproducible:

```ini
target_prep_dir = /data/XTT/01.prepare
```

and the reference:

```ini
reference_gff = /data/Eru.gff3
reference_cds = /data/Eru.cds
reference_pep =            ; optional; translated from reference_cds if omitted
```

### `[chromosomes]`
A polyploid is described as a set of **homology groups**. Two ways to say so:

**standard** — chromosomes are named `<prefix><group><suffix>`, e.g. `Chr1A … Chr10L`:

```ini
mode        = standard
prefix      = Chr
group_start = 1
group_end   = 10
subgenomes  = A,B,C,D,E,F,G,H,I,J,K,L
```

**group_file** — anything else; list the groups explicitly
(see [`groups.example.tsv`](groups.example.tsv)):

```ini
mode       = group_file
group_file = ./groups.tsv
```

> Groups do **not** have to contain the same number of chromosomes. Sugarcane-like
> auto-allopolyploids regularly have 5–12 homoeologs per group, and the pipeline
> handles that natively.

### `[reference]`
Which reference chromosome defines the frame for each group:

```ini
group_ref_map = Group_01:Chr01, Group_02:Chr02, ...
```

Leave it empty to match by number (group 1 → `Chr01`, …).

### `[output]`
```ini
work_dir = ./polyallele_work      ; scratch
out_dir  = ./polyallele_out       ; final tables
prefix   = polyallele
```

### `[runtime]`
```ini
conda_env  = polyallele
scheduler  = none | pbs | slurm
pbs_ppn    = 16
pbs_walltime = 12:00:00
```

`none` runs in the current shell (a laptop, or a node you are already on).
`pbs` / `slurm` writes a job script into `work_dir` and submits it.

---

## How the algorithm works

```
   target chromosomes                 reference (diploid)
          │                                   │
          │   (1) PREP                        │  split into per-chromosome
          │       split / extract             │  .cds  .bed  .pep
          │                                   │
          ▼                                   ▼
   ┌──────────────────────────────────────────────────────┐
   │ (2) PROJECT   JCVI MCScan, target chromosome vs its  │
   │               own reference chromosome, pairwise     │
   └──────────────────────────────────────────────────────┘
                          │  .anchors  (gene ↔ gene, positional)
                          ▼
   ┌──────────────────────────────────────────────────────┐
   │ (3) INTERP    walk each chromosome in gene order.    │
   │               Two consecutive anchors bracket an     │
   │               interval; the genes between them are   │
   │               assigned by ORDER-PRESERVING          │
   │               interpolation onto the reference loci  │
   │               in that interval.                      │
   │                                                      │
   │               · inversions handled (reverse walk)    │
   │               · no room in the reference -> a        │
   │                 SYNTHETIC interval slot, so no cell  │
   │                 ever receives two genes              │
   └──────────────────────────────────────────────────────┘
                          │
                          ▼
   ┌──────────────────────────────────────────────────────┐
   │ (4) REFINE    a gene alone in a small row is pulled  │
   │               into a nearby rich row -- but ONLY if  │
   │               (a) that row is within `window`        │
   │                   reference loci on the SAME         │
   │                   reference chromosome, AND          │
   │               (b) BLASTP supports it.                │
   │                                                      │
   │               Similarity never widens the window;    │
   │               it only decides inside it.             │
   └──────────────────────────────────────────────────────┘
                          │
                          ▼
   ┌──────────────────────────────────────────────────────┐
   │ (5) REPORT    wide matrix + long table + unplaced    │
   └──────────────────────────────────────────────────────┘
```

The 1:1 property — *at most one gene per (locus, chromosome) cell* — is
**guaranteed by construction** in stage 3: a candidate locus already taken by this
chromosome is never overwritten; the gene is pushed to its own synthetic slot.
The only multi-gene cells that survive are those where two **anchored** genes
project onto the same reference locus, which is a genuine copy-number signal.

---

## Output files

| file | contents |
|---|---|
| `<prefix>.allele_matrix.tsv` | **the deliverable.** one row per locus; `Allele_*` columns; `Extra_Copies`; `Extra_Classes`; `NCopies`; `Synthetic`; `Evidence`; `Merged` |
| `<prefix>.allele_long.tsv` | one line per gene: `LocusID, HomologyGroup, Ref_Gene, Ref_Locus, Subgenome, Chromosome, GeneID, Rank, CopyClass, Source, Interval, Synthetic` |
| `<prefix>.synthetic_rows.tsv` | the reference-bracketed rows (no `Ref_Gene`) in long form — positional annotations, *not* verified loci |
| `<prefix>.unplaced_genes.tsv` | genes that entered no row at all |
| `<prefix>.unplaced_matrix.tsv` | **those same genes, arranged as a matrix** — rows are `>= identity` components computed over the unplaced genes only, columns are chromosomes. A gene with no partner gets its own row, and *that* row really is a single-copy gene |
| `<prefix>.locus2gene.tsv` | `LocusID <TAB> GeneID` plus chromosome, column, `Role` and `Source` — table genes only |
| `<prefix>.gene2copyclass.tsv` | **every** gene, labelled `single_copy` or `multi_copy` from sequence similarity alone |
| `<prefix>.stats.txt` | summary counts, including the copy-identity tally |
| `polyallele.log` | full run log |

### `Evidence` — how strong is the positional claim?

A row states *which ancestral locus* it is. How much that statement is worth
depends entirely on how the locus was reached, so it is a column:

| `Evidence` | meaning | genes in the reference run |
|---|---|---|
| `anchor` | the row holds at least one direct collinear anchor | 187,902 |
| `interval` | no anchor, but the locus is a real reference gene, reached by interpolation inside a bracket | 82,012 |
| `slot` | the row is a synthetic slot: the reference interval is known, the locus is not | 42,725 |

`Merged` counts how many genes in the row were brought in by the `merge` stage.

### The two copy tables — and what they are *not*

`locus2gene.tsv` groups by **ancestral locus**. `gene2copyclass.tsv` groups by
**sequence similarity alone**, with position deliberately ignored. They answer
different questions and must not be swapped:

| | `locus2gene.tsv` | `gene2copyclass.tsv` |
|---|---|---|
| a cluster is | one ancestral locus | one similarity component |
| answers | who is allelic to whom | which genes have close relatives |
| for allele-level quantification | ✅ | ❌ — one component can span several loci |

A component in `gene2copyclass.tsv` is a **transitive closure**: gene A joins B
and B joins C at `>= min_identity`, so A and C share a cluster even if they are
only 70% identical.

**Build it from BLASTP, not from DIAMOND `--sensitive`.** That was measured:
against BLASTP on the same data, **24.58% of DIAMOND's `single_copy` calls were
wrong** — the gene did have a `>=90%` relative. The whole point of the column is
to say "nothing else in this genome is like this gene", so sensitivity is the
thing being measured.

Two further cautions from the reference run:

* **"Within group" is not "whole genome".** Of 42,659 genes that were alone inside
  their homology group, **24,045 (56%) had a `>=90%` relative in a DIFFERENT
  group.** Screening only within groups would have called all 42,659 single copy.
* **A cross-group relative is not proof of an error in the grouping.** That
  finding is large enough that "different homology group" cannot be treated as
  automatically impossible when adjudicating a merge; it is strong evidence, not
  proof.

The reference run ends with **18,614 `single_copy`** and **317,597 `multi_copy`
genes** out of 336,211.

### `unplaced_matrix.tsv`

The genes that entered no row are not a shapeless list, and they are not all
single copy either. Arranged the same way as the main matrix, with rows formed
from `>= identity` components over the unplaced genes **only**:

| `Evidence` | meaning | reference run |
|---|---|---|
| `no_partner` | the component has one member — nothing else among the unplaced is its allele | 16,026 rows |
| `similarity_pair` | two or more unplaced genes are each other's alleles, so neither is single copy even though neither is in the table | 2,759 rows |

`Chromosomes` lists the chromosomes the row spans; `NCopies` the genes in it.

`Source` in the long table tells you how each gene got its locus:

| `Source` | meaning | positional confidence |
|---|---|---|
| `anchor` | direct collinear anchor to the reference | highest |
| `copy_at_anchored_locus` | both flanking anchors agree on locus L and this gene lies physically between them → an additional copy **at L** | high |
| `interp_exact` | interval had exactly as many target genes as reference loci → forced 1:1 | high |
| `interp_proportional` | fewer target genes than reference loci → order alone cannot decide which | medium — check `Ref_Gene` |
| `interp_similarity` | as above, resolved with BLASTP inside the interval | medium-high |
| `interp_synthetic` | no room in the reference → own interval slot | interval is reliable, locus is not |
| `merged_homology` | satellite pulled back into a nearby rich row | depends on the target row |

### Copy identity — what kind of paralog?

An extra copy is not just "an extra". Its **physical distance from the cell's
representative** says which duplication mechanism produced it:

| `CopyClass` | distance | usually means |
|---|---|---|
| `tandem_copy` | ≤ `output.tandem_kb` (default 10 kb) | tandem array — unequal crossing over / slippage |
| `proximal_copy` | ≤ `output.proximal_kb` (default 1000 kb) | local / segmental duplication |
| `dispersed_copy` | larger | dispersed duplication, or a transposed copy |
| `representative` | rank 1 in the cell | the copy that stands for the cell |

That is the difference between reporting "copy number 3" and reporting
"one homoeolog plus one tandem array plus one dispersed segmental duplicate".

> **Scope of the claim.** Distance classifies *the physical relationship*, which is
> solid. Inferring *mechanism* (tandem vs segmental vs whole-genome duplication) from
> distance alone is an interpretation — a syntenic block-level check is needed to
> call segmental duplication with confidence. Use `CopyClass` as the physical fact,
> and confirm the mechanism before writing it into a paper.

---

## How to read the output

**Downstream quantification — the intended use.**

```
one row            = one ancestral locus          -> never merge rows
one cell (row,col) = one subgenome's copy(ies)    -> merge inside a cell
columns            = subgenomes                   -> NEVER merge across columns
```

That is the collapse map the pipeline exists to deliver:

1. **Within a cell**, sum the representative and its `Extra_Copies`. They are
   copies of the same locus on the same chromosome — indistinguishable to a
   quantifier and equivalent biologically.
2. **Across columns**, do not collapse. Homoeolog ambiguity is precisely the
   signal you want; let the quantifier's EM resolve it.
3. **Across rows**, never collapse. Different rows are different genes.

**Filtering by confidence.** For the cleanest analysis use `Source = anchor` plus
`interp_exact` plus `interp_similarity`. Keep `interp_proportional` and synthetic
rows for coverage-sensitive work, and say so in your methods.

---

## Quality control

`bin/qc_matrix.py` reports the numbers that decide whether a matrix is usable:

```
[ 2 ] THE 1:1 PROPERTY
  cells with exactly 1 gene      : ...
  inflation (genes / cells)      : ...
[ 3 ] SUBGENOME OCCUPANCY          <- should peak near the group's chromosome count
[ 5 ] OVER-MERGE CHECK             <- reference loci spread over >1 row: must be 0
[ 6 ] GENE UNIQUENESS              <- genes appearing twice: must be 0
[ 7 ] SPAN PURITY                  <- are multi-gene cells tandem or dispersed?
```

**Validate against a curated gene family.** The single most informative test is to
take a family you have already annotated by hand (a dozen members is plenty) and
check two things:

- **over-merge** — does any row contain two members of *different* families? (must be 0)
- **over-split** — how many rows does one family member occupy, counting only the
  homology group it actually lives in? (ideal 1.0)

On the sugarcane dataset this pipeline gave over-merge = 0 and a within-group
over-split of 1.42, against 2 and 1.25 respectively for a similarity-clustering
baseline. The trade is deliberate: the pipeline buys the fatal error away at the
cost of a little fragmentation.

---

## Parameters that matter

| parameter | section | effect |
|---|---|---|
| `cscore` | `[projection]` | lower → more anchors → higher coverage, noisier. `0.7` for a congeneric reference, `0.9–0.99` for a distant one |
| `anchors_dir` | `[input]` | reuse already-computed `.anchors` files instead of recomputing the projection |
| `max_interval` | `[interpolation]` | intervals wider than this many reference loci are left unplaced instead of guessed |
| `max_copies_at_anchored_locus` | `[interpolation]` | cap on how many genes one anchored locus may absorb before the surplus is pushed to its own slot |
| `similarity_first` | `[refinement]` | resolves proportionally-placed genes with BLASTP *inside* the interval |
| `rich_min` / `sparse_max` | `[refinement]` | which rows may absorb / be absorbed |
| `window` | `[refinement]` | how far, in reference loci, a merge is allowed to reach. **The main safety dial** |
| `min_identity` / `min_alignment` | `[refinement]` | BLASTP stringency for a merge |
| `tandem_kb` / `proximal_kb` | `[output]` | where the boundaries between `tandem_copy` / `proximal_copy` / `dispersed_copy` are drawn |

Tuning `cscore` upward lowers coverage but raises confidence; if coverage becomes
the binding constraint, prefer a **closer reference genome** over loosening
`cscore` — it is by far the larger lever.

---

## Limitations — read this before publishing

1. **Coverage is not 100%, and it should not be reported as if it were.**
   Genes outside any collinear block cannot be placed by position, and inventing a
   placement for them is what makes similarity pipelines over-merge. On the
   sugarcane test ~92% of genes were placed, ~68% into reference-anchored rows.

2. **Synthetic rows are positional annotations, not verified loci.** They carry no
   `Ref_Gene`. Their *interval* is reliable (roughly 85% of the time the true locus
   lies inside it); their *exact* locus is not.

3. **`interp_proportional` is the weakest placement class.** Order alone cannot say
   which of several reference loci a gene belongs to; expect the exact locus to be
   wrong in a substantial fraction of cases, though almost always within a couple of
   loci of the truth. Enable `similarity_first` to improve it.

4. **Coverage is not uniform across chromosomes.** Chromosomes carrying large
   lineage-specific insertions, or poor assemblies, anchor badly. `qc_matrix.py`
   and the per-chromosome lines in the log will show this; investigate before
   trusting those groups.

5. **The pipeline needs a good diploid reference.** If no clean diploid relative
   exists for your system, this approach — not just this implementation — is the
   wrong tool.

6. **Rows with few copies are weak groups.** A row holding 2 of a possible 12
   homoeologs is a legitimate finding (gene loss) but a poor quantification unit.
   Consider requiring a minimum occupancy downstream.

7. **This is a hypothesis generator.** Validate with a curated family before
   drawing biological conclusions.

---

## Resuming after an interruption

Every stage is idempotent:

```bash
bash run_polyallele.sh config.ini project    # skips anchor files that already exist
bash run_polyallele.sh config.ini refine     # re-runs the cheap upstream stages
bash run_polyallele.sh config.ini report     # rebuilds the tables
```

Running a later stage automatically runs its prerequisites. The expensive stage
(`project`, the JCVI pairwise alignments) skips any pair whose `.anchors` file is
already present, so it is safe to re-run after a walltime kill.

---

## Citation

> Manuscript in preparation. Please cite this repository and the reference genome
> you used until a paper is available.

If you use the JCVI/MCScan implementation, cite:

> Tang H. et al. *Synteny and collinearity in plant genomes.* Science 320, 486–488 (2008).

---

## License

Non-commercial research licence — see [LICENSE](LICENSE).
Academic and research use is free; commercial use is prohibited.
