# Dundee cluster interactive instructions (no shared filesystem)

These notes describe how to run the **myRna-seq** Snakemake workflow on the Dundee
cluster **inside a single interactive job**. Because the cluster has **no shared
filesystem**, you must stage the project, references, and FASTQs on the compute
node and write outputs to node-local storage. At the end, sync the results back
to your persistent storage.

## 1) Start an interactive session

Follow the same interactive workflow used for other Dundee pipelines, then
request resources with the SGE flags used by this project.

```bash
# 1. Connect to the cluster
ssh compute.dundee.ac.uk

# 2. Start a screen session (protects against connection loss)
screen -S snakemake_job

# If connection drops, reconnect with:
# ssh compute.dundee.ac.uk
# screen -r snakemake_job

# 3. Request an interactive session with the SGE settings used by this pipeline
qrsh -pe smp 40 \
  -adds l_hard local_free 200G \
  -mods l_hard m_mem_free 20G \
  -adds l_hard avx 1
```

Once the session starts you will be on a compute node.

## 2) Stage the project and inputs on local scratch

Navigate to your lab folder, clone the repository (first time only), and set up
your environment (first time only). Then copy data to node-local scratch for
the run.

```bash
# 4. Navigate to your lab folder
cd /cluster/majf_lab/mtinti  # Replace with your lab folder path

# 5. Clone the repository (this step is done only once)
git clone https://github.com/mtinti/myRna-seq.git

# 6. Create and activate conda environment (this step is done only once)
conda create -n snakemake snakemake
# 6b. (This step is every time)
conda activate snakemake

# 7. Navigate to the repository
cd myRna-seq
# 7b. git pull (this step is every time)
```

Pick a run directory on the node (`$TMPDIR` is node-local) and copy the
pipeline repository plus your input data there.

```bash
RUN_DIR="$TMPDIR/myRna-seq-$USER-$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RUN_DIR"

# Copy the pipeline repository (exclude git metadata to keep it small)
rsync -a --exclude '.git' /cluster/majf_lab/mtinti/myRna-seq/ "$RUN_DIR/myRna-seq/"

# Copy input FASTQs and references (examples)
rsync -a /path/to/fastqs/ "$RUN_DIR/fastqs/"
rsync -a /path/to/reference/ "$RUN_DIR/reference/"
```

> **Tip:** If your login node storage is not visible from the compute node,
> replace the `rsync` steps with whatever transfer method you use on Dundee
> (e.g., `scp` from a data transfer node).

## 3) Configure the workflow for node-local paths

Set the paths in `config.yaml` to match where you staged the data. You can either
edit the file directly or override values at runtime.

Example override at runtime:

```bash
cd "$RUN_DIR/myRna-seq"

snakemake --cores 40 \
  --config \
  processing_dir="$RUN_DIR/processing" \
  results_dir="$RUN_DIR/results" \
  benchmark_dir="$RUN_DIR/benchmarks" \
  reference_fasta="$RUN_DIR/reference/genome.fa" \
  gtf_file="$RUN_DIR/reference/annotation.gtf" \
  samples_csv="$RUN_DIR/samples.csv"
```

## 4) Sync results back to persistent storage

When the workflow finishes, copy outputs back to your long-term storage:

```bash
rsync -a "$RUN_DIR/results/" /path/to/output/results/
rsync -a "$RUN_DIR/benchmarks/" /path/to/output/benchmarks/
```

You can also archive logs if needed:

```bash
rsync -a "$RUN_DIR/.snakemake/log/" /path/to/output/snakemake-logs/
```

---

### Notes

- This approach runs the full workflow on **one node**; avoid Snakemake cluster
  submission (`--cluster`) because compute nodes do not share a filesystem.
- Set `max_cores` and the `cores_*` values in `config.yaml` to match your job
  allocation if you override these defaults.
