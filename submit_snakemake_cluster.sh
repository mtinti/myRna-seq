#!/bin/bash
#$ -adds l_hard local_free 200G
#$ -mods l_hard m_mem_free 20G
#$ -adds l_hard avx 1
#$ -cwd
#$ -V
#$ -j y
#$ -N snakemake_rnaseq
#$ -o snakemake_rnaseq_$JOB_ID.log
#$ -pe smp 40

# Exit on error and undefined variables
set -e
set -u

# ── Activate conda ──────────────────────────────────────────────
eval "$(conda shell.bash hook)"
conda activate snakemake

# ── User-configurable variables ─────────────────────────────────
CORES="${CORES:-40}"
CONFIGFILE="${CONFIGFILE:-config.yaml}"

# ── Override processing_dir to high-performance scratch ─────────
# The pipeline copies final outputs from processing_dir to results_dir
# internally, so results land directly in results/ (persistent storage).
# No rsync needed — this is the same approach as the interactive workflow.
PROCESSING_DIR="${TMPDIR:-/tmp}/processing"

# ── Run Snakemake ───────────────────────────────────────────────
snakemake \
  --cores "$CORES" \
  --configfile "$CONFIGFILE" \
  --use-conda \
  --config "processing_dir=$PROCESSING_DIR"
