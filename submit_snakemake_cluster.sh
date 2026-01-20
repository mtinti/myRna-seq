#!/bin/bash
#$ -adds l_hard local_free 200G
#$ -mods l_hard m_mem_free 20G
#$ -adds l_hard avx 1
#$ -cwd
#$ -V
#$ -j y
#$ -N snakemake_oligo
#$ -o snakemake_oligo_errors_$JOB_ID
#$ -pe smp 40

# Exit on error and undefined variables
set -e
set -u

WORKDIR="${SGE_O_WORKDIR:-$PWD}"
JOB_ID="${JOB_ID:-manual}"
RUN_ROOT="${TMPDIR:-/tmp}"
RUN_DIR="$RUN_ROOT/myRna-seq-${USER}-${JOB_ID}"
RESULTS_DEST="${RESULTS_DEST:-$WORKDIR}"

CORES="${CORES:-40}"
CONFIG_OVERRIDES=(
  "processing_dir=$RUN_DIR/processing"
  "results_dir=$RUN_DIR/results"
  "benchmark_dir=$RUN_DIR/benchmarks"
)

mkdir -p "$RUN_DIR"

rsync -a \
  --exclude '.git' \
  --exclude 'processing' \
  --exclude 'results' \
  --exclude 'benchmarks' \
  "$WORKDIR/" "$RUN_DIR/myRna-seq/"

cd "$RUN_DIR/myRna-seq"

SNK_ARGS=(--cores "$CORES" --config "${CONFIG_OVERRIDES[@]}")

snakemake "${SNK_ARGS[@]}"

mkdir -p "$RESULTS_DEST/results" "$RESULTS_DEST/benchmarks" "$RESULTS_DEST/snakemake-logs"

rsync -a "$RUN_DIR/results/" "$RESULTS_DEST/results/"
rsync -a "$RUN_DIR/benchmarks/" "$RESULTS_DEST/benchmarks/"

if [[ -d "$RUN_DIR/.snakemake/log" ]]; then
  rsync -a "$RUN_DIR/.snakemake/log/" "$RESULTS_DEST/snakemake-logs/"
fi

echo "Run directory: $RUN_DIR"
echo "Results synced to: $RESULTS_DEST"
