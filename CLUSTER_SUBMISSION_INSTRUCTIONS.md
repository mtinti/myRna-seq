# Dundee cluster batch submission instructions (no shared filesystem)

These instructions submit **myRna-seq** as a single batch job using `qsub` and a
wrapper script (`submit_snakemake_cluster.sh`). The script runs directly from
your persistent lab folder — it only redirects `processing_dir` to the
high-performance scratch drive (`$TMPDIR`) at execution time. The pipeline
copies final results back to `results/` automatically, so no manual rsync is
needed.

The wrapper automatically handles:
- Activating the conda environment
- Overriding `processing_dir` to `$TMPDIR/processing` (any value in the config file is replaced)
- Running Snakemake with `--use-conda` and your config file

## 1) Set up the environment (one-time)

```bash
# 1. Connect to the cluster
ssh compute.dundee.ac.uk

# 2. Navigate to your lab folder
cd /cluster/majf_lab/mtinti  # Replace with your lab folder path

# 3. Clone the repository (this step is done only once)
git clone https://github.com/mtinti/myRna-seq.git

# 4. Create the conda environment (this step is done only once)
conda create -n snakemake snakemake
```

## 2) Prepare the input data

```bash
# 5. Navigate to the repository
cd myRna-seq

# 6. Pull the latest changes (every time)
git pull

# 7. Create the input data directory
mkdir -p indata
```

Copy your FASTQ files into a sample-specific subfolder. For example:

```
indata/TbRIT-8166/V350168884_L02_B5GTBRjalrRAACA-4_1.fq.gz
indata/TbRIT-8166/V350168884_L02_B5GTBRjalrRAACA-4_2.fq.gz
indata/TbRIT-8166/V350194510_L02_B5GTBRjalrRAACA-4_1.fq.gz
indata/TbRIT-8166/V350194510_L02_B5GTBRjalrRAACA-4_2.fq.gz
```

## 3) Prepare the sample CSV file

Create a sample CSV (e.g. `samples_rit.csv`). It should look like:

```csv
sample_name,read_type,source_type,file_path_1,file_path_2,checksum_1,checksum_2,run_tag
SAMPLE_rit1,paired,local,indata/TbRIT-8166/V350168884_L02_B5GTBRjalrRAACA-4_1.fq.gz,indata/TbRIT-8166/V350168884_L02_B5GTBRjalrRAACA-4_2.fq.gz,db36979b3164a2a29c9af85e6b3072fd,01aad38e21bc8ea452612b10037a58a7,TbRIT-8166
SAMPLE_rit2,paired,local,indata/TbRIT-8166/V350194510_L02_B5GTBRjalrRAACA-4_1.fq.gz,indata/TbRIT-8166/V350194510_L02_B5GTBRjalrRAACA-4_2.fq.gz,876109f629db0fef6c93918137d467a0,6871417faaf83beae9c72051a12519ce,TbRIT-8166
```

Key notes about the sample CSV:

- **`sample_name`** must be unique for each row (e.g. `SAMPLE_rit1`, `SAMPLE_rit2`).
- **`run_tag`** groups samples that belong to the same experiment; samples sharing a `run_tag` will be merged together.
- **`checksum_1` / `checksum_2`** should contain the md5sum from BIG, to verify that the files you are analysing match what BIG produced.
- Please read the myRna-seq README for a full description of all available columns.

## 4) Create and edit the config file

```bash
# 8. Make a copy of the default config
cp config.yaml config_rit.yaml
```

This copy serves both as the configuration for your run and as a record of how the pipeline was executed. Parameters can also be overridden at the command line (see the myRna-seq README for details).

Edit `config_rit.yaml` and apply the following changes:

### a) Processing directory (handled automatically)

You do **not** need to change `processing_dir` in the config file. Whatever value is set in the config is automatically replaced at execution time with the compute node's high-performance scratch drive (`$TMPDIR/processing`). The pipeline processes data there and copies results back to `results/` on persistent storage by itself.

### b) Set the reference genome and annotation

Replace the default reference paths:

```yaml
reference_fasta: "reference/genome_427/tb427.fa"
gtf_file: "reference/genome_427/tb427.gtf"
```

with the assembly you want to use. For example, the TriTrypDB v68 TREU927 assembly:

```yaml
reference_fasta: "reference/genome_927_68/TriTrypDB-68_TbruceiTREU927_Genome.fa"
gtf_file: "reference/genome_927_68/TriTrypDB-68_TbruceiTREU927.gff"
```

Notes:
- The `reference/genome_927_68/` directory should sit at the same level as `config_rit.yaml`.
- The FASTA file **must** end with `.fa` — TriTrypDB files may need to be renamed from `.fasta`.
- The pipeline can automatically generate a GTF file from a GFF, so providing a `.gff` is fine.

### c) Point to your sample CSV

Change:

```yaml
samples_csv: "samples_full.csv"
```

to:

```yaml
samples_csv: "samples_rit.csv"
```

### d) Enable BAM copying

Change:

```yaml
copy_bam: False
```

to:

```yaml
copy_bam: True
```

This copies the BAM files from the processing directory back into the results folder. The BAM files will serve as input for further processing with **myBarcode-Seq**.

## 5) Submit the job

```bash
# 9. Activate the conda environment
conda activate snakemake

# 10. Submit the job with your custom config file
CONFIGFILE=config_rit.yaml qsub submit_snakemake_cluster.sh
```

To use the default `config.yaml` instead, simply run:

```bash
qsub submit_snakemake_cluster.sh
```

Check the job queue and review the log once it finishes:

```bash
# Check job status
qstat

# Review the log (replace JOBID with your job ID)
cat snakemake_rnaseq_JOBID.log
```

## 6) Results

When the pipeline finishes, results are already in `results/` inside your project directory — no manual copying needed. The pipeline handles moving outputs from scratch back to persistent storage internally.

Contents of `results/`:
- Feature counts and RPKM values
- Coverage tracks (bigWig)
- QC reports
- BAM files (if `copy_bam: True`)
- Benchmark summaries

---

### Notes

- Because the cluster has **no shared filesystem**, the workflow runs on a
  **single node**. Do not use Snakemake's cluster submission mode.
- The wrapper script uses `--use-conda` so that rule-level conda environments
  are created and activated automatically.
- The `processing_dir` value in your config file is always overridden by the
  script — you can ignore it or leave it at the default.
- `CONFIGFILE=config_rit.yaml qsub ...` sets an environment variable that the
  script reads via `${CONFIGFILE:-config.yaml}`. The `#$ -V` directive in the
  script exports all environment variables to the job, so the value you set
  before `qsub` is available inside the running job.
- To change resource requests (cores, memory, disk), edit the SGE directives
  (`#$` lines) at the top of `submit_snakemake_cluster.sh`. Make sure `CORES`
  matches `#$ -pe smp`.

> **Tip:** Keep separate config and sample files for each experiment. This way the pipeline run is fully reproducible without needing to store intermediate files.