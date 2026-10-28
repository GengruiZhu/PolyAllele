#!/usr/bin/env bash
# =============================================================================
#  PolyAllele v4 -- one command to run everything.
#
#      bash run_polyallele.sh              # uses ./config.ini
#      bash run_polyallele.sh my.ini       # uses another config
#      bash run_polyallele.sh config.ini project    # a single stage
#
#  Everything is driven by config.ini: inputs, parameters, outputs and how to
#  launch (directly, or through PBS / SLURM).
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:-$HERE/config.ini}"
STAGE="${2:-all}"

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: config file not found: $CONFIG" >&2
    exit 2
fi
CONFIG="$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")"

# --- pick a python ----------------------------------------------------------
PY="${PYTHON:-}"
if [ -z "$PY" ]; then
    for c in python3 python; do
        if command -v "$c" >/dev/null 2>&1; then PY="$c"; break; fi
    done
fi
[ -z "$PY" ] && { echo "ERROR: no python3 on PATH" >&2; exit 3; }

# --- read the runtime section out of config.ini ------------------------------
eval "$("$PY" "$HERE/bin/polyallele.py" --config "$CONFIG" --emit-shell)"

# --- activate the conda environment, if one was named ------------------------
# conda's activate.d scripts reference variables that may be unset, which trips
# `set -u`.  Relax nounset for the activation only.
if [ -n "$POLYALLELE_CONDA_ENV" ]; then
    set +u
    if command -v conda >/dev/null 2>&1; then
        eval "$(conda shell.bash hook)"
    elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        . "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
        . "$HOME/anaconda3/etc/profile.d/conda.sh"
    elif [ -x "$HOME/miniconda3/bin/conda" ]; then
        eval "$("$HOME"/miniconda3/bin/conda shell.bash hook)"
    fi
    if command -v conda >/dev/null 2>&1; then
        conda activate "$POLYALLELE_CONDA_ENV" 2>/dev/null || \
            echo "WARNING: could not activate conda env '$POLYALLELE_CONDA_ENV'; using $PY"
    else
        echo "WARNING: conda not found; using $PY"
    fi
    set -u
fi

echo "============================================================"
echo " PolyAllele v4"
echo "   config    : $CONFIG"
echo "   stage     : $STAGE"
echo "   python    : $(command -v "$PY")"
echo "   scheduler : $POLYALLELE_SCHEDULER"
echo "============================================================"

# --- direct run --------------------------------------------------------------
if [ "$POLYALLELE_SCHEDULER" = "none" ]; then
    "$PY" "$HERE/bin/polyallele.py" --config "$CONFIG" "$STAGE"
    STATUS=$?
    [ "$STATUS" -ne 0 ] && exit "$STATUS"

    # --- the tables that sit next to the matrix ------------------------------
    # Kept in a separate script on purpose: they only READ what the pipeline
    # wrote, so a problem here can never damage the matrix, and they can be
    # re-run on their own after changing a threshold.
    if [ "${POLYALLELE_EXTRAS:-0}" = "1" ] && \
       { [ "$STAGE" = "all" ] || [ "$STAGE" = "report" ]; }; then
        if [ -z "${POLYALLELE_ALLVALL:-}" ]; then
            echo "NOTE: copyclass.allvall is empty -- the extra tables need a"
            echo "      genome-wide all-vs-all screen; skipping them."
        elif [ ! -s "$POLYALLELE_ALLVALL" ]; then
            echo "NOTE: $POLYALLELE_ALLVALL is missing or empty; skipping the extra tables."
        else
            "$PY" "$HERE/bin/extra_tables.py" \
                --matrix   "$POLYALLELE_OD/$POLYALLELE_PREFIX.allele_matrix.tsv" \
                --long     "$POLYALLELE_OD/$POLYALLELE_PREFIX.allele_long.tsv" \
                --unplaced "$POLYALLELE_OD/$POLYALLELE_PREFIX.unplaced_genes.tsv" \
                --allvall  "$POLYALLELE_ALLVALL" \
                --out-dir  "$POLYALLELE_OD" \
                --prefix   "$POLYALLELE_PREFIX" \
                --identity "$POLYALLELE_IDENTITY" \
                --competition-identity "$POLYALLELE_COMPID" \
                --min-alignment "$POLYALLELE_MINLEN" || \
                echo "WARNING: extra_tables.py failed; the matrix itself is unaffected"
        fi
    fi
    exit 0
fi

# --- scheduled run -----------------------------------------------------------
mkdir -p "$POLYALLELE_WORKDIR"
JOB="$POLYALLELE_WORKDIR/run_polyallele.job.sh"
NPROC="$POLYALLELE_PPN"

if [ "$POLYALLELE_SCHEDULER" = "pbs" ]; then
    cat > "$JOB" <<EOF
#!/bin/bash
#PBS -N $POLYALLELE_JOBNAME
#PBS -l nodes=1:ppn=$NPROC,walltime=$POLYALLELE_WALLTIME
#PBS -q $POLYALLELE_QUEUE
#PBS -o $POLYALLELE_WORKDIR/polyallele.job.log
#PBS -j oe
cd "\$PBS_O_WORKDIR" 2>/dev/null || cd "$HERE"
bash "$HERE/run_polyallele.sh" "$CONFIG" "$STAGE"
EOF
    echo ">>> qsub $JOB"
    qsub "$JOB"
elif [ "$POLYALLELE_SCHEDULER" = "slurm" ]; then
    cat > "$JOB" <<EOF
#!/bin/bash
#SBATCH --job-name=$POLYALLELE_JOBNAME
#SBATCH --cpus-per-task=$NPROC
#SBATCH --time=$POLYALLELE_WALLTIME
#SBATCH --output=$POLYALLELE_WORKDIR/polyallele.job.log
cd "\$SLURM_SUBMIT_DIR" 2>/dev/null || cd "$HERE"
bash "$HERE/run_polyallele.sh" "$CONFIG" "$STAGE"
EOF
    echo ">>> sbatch $JOB"
    sbatch "$JOB"
else
    echo "ERROR: unknown scheduler '$POLYALLELE_SCHEDULER' (use none / pbs / slurm)" >&2
    exit 4
fi
