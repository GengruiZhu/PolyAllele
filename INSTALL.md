# Installation

PolyAllele v4 is pure Python plus four external command-line tools. There is no
build step and nothing to compile.

---

## 1. Requirements

| requirement | why |
|---|---|
| Python ≥ 3.8 | the pipeline itself (standard library only) |
| `numpy` | `bin/qc_matrix.py` |
| `jcvi` (Python package) | MCScan collinear projection — `jcvi.compara.catalog` |
| `last` (`lastal`, `lastdb`) | the aligner JCVI uses for that projection |
| `blast` (`makeblastdb`, `blastp`) | the refinement stage |
| `gffread` | only when `input.target_prep_dir` is **not** supplied |

Disk: the working directory holds per-chromosome sequences, one `.last` file per
chromosome pair and BLAST databases. Budget roughly 3–5× the size of your genome
FASTA. For a ~1.5 Gb genome that is ~5–8 GB.

Memory: under 4 GB for genomes up to a few Gb; the BLAST database of a large
reference protein set dominates.

---

## 2. Conda (recommended)

```bash
conda env create -f environment.yml
conda activate polyallele
```

Then verify:

```bash
python      --version
lastal      --version | head -1
blastp      -version | head -1
gffread     --version
python -c "import jcvi, numpy; print('python deps OK')"
python bin/polyallele.py --help
```

`environment.yml` already requests everything. On a machine without internet
access, ask your administrator to pre-populate the environment, or create it from
a local channel mirror — the pipeline itself never downloads anything.

---

## 3. Manual installation (no conda)

Install the four tools by whatever means your system provides, then make sure they
are on `PATH`:

```bash
for t in lastal makeblastdb blastp gffread; do command -v $t || echo "MISSING: $t"; done
python -c "import jcvi, numpy" || pip install jcvi numpy
```

`jcvi` needs no configuration file for this pipeline. If your `jcvi` build rejects
an option, note that PolyAllele probes for `--cpus` and adapts automatically.

---

## 4. Verifying the installation on a small dataset

Run the `prep` stage only, then look at the log:

```bash
bash run_polyallele.sh config.ini prep
```

A healthy `prep` prints, at the end:

```
[prep] N homology groups, M target chromosomes
[prep]   Group_01 -> reference Chr01: 11 chromosomes Chr1A..Chr1K
...
[prep] manifest: .../polyallele_work/manifest.tsv
```

Common failures:

| message | cause |
|---|---|
| `reference_gff and reference_cds are required` | `[input]` not filled in |
| `gffread not found in PATH` | install gffread, or supply `target_prep_dir` |
| `no CDS for target chromosome X` | the chromosome has no `.cds`, or its name does not match the group definition |
| `cannot determine the reference chromosome for Group_XX` | the group name carries no number — fill in `reference.group_ref_map` |

---

## 5. Running on a cluster

Set the scheduler block in `config.ini`:

```ini
[runtime]
conda_env    = polyallele
scheduler    = pbs          ; or slurm
pbs_queue    = comput
pbs_ppn      = 16
pbs_walltime = 12:00:00
```

`bash run_polyallele.sh` now writes a job script into `work_dir` and submits it.
Nothing runs on the login node.

The pipeline is **resumable**: if the walltime kills it, re-run the same command and
the `project` stage will skip every chromosome pair whose `.anchors` file already
exists. For large genomes, submit `project` on its own first, then the rest:

```bash
bash run_polyallele.sh config.ini project
bash run_polyallele.sh config.ini report
```

**Cluster note — nothing runs on the login node.** If you are wiring this into an
agent or a batch harness, keep that rule: `run_polyallele.sh` with
`scheduler = none` executes *wherever you call it*, so only use `none` on a node
you are entitled to compute on.

---

## 6. Uninstalling

Delete the directory tree:

```bash
rm -rf polyallele_work polyallele_out
```

Nothing is installed outside the working directory and the conda environment.
