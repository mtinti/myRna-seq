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

## 3) Prepare the input data

Create an input data directory and copy your FASTQ files into it:

```bash
# 8. Create the input data directory
mkdir -p indata
```

Copy your FASTQ files into a sample-specific subfolder. For example:

```
indata/TbRIT-8166/V350168884_L02_B5GTBRjalrRAACA-4_1.fq.gz
indata/TbRIT-8166/V350168884_L02_B5GTBRjalrRAACA-4_2.fq.gz
indata/TbRIT-8166/V350194510_L02_B5GTBRjalrRAACA-4_1.fq.gz
indata/TbRIT-8166/V350194510_L02_B5GTBRjalrRAACA-4_2.fq.gz
```

## 4) Identify the high-performance scratch directory

```bash
# 9. Find out where the high-performance local drive is on your node
echo $TMPDIR/
# You should see something like: /tmp/3544738.1.rhel9.q/
```

## 5) Prepare the sample CSV file

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

## 6) Create and edit the config file

```bash
# 10. Make a copy of the default config
cp config.yaml config_rit.yaml
```

This copy serves both as the configuration for your run and as a record of how the pipeline was executed. Parameters can also be overridden at the command line (see the myRna-seq README for details).

Edit `config_rit.yaml` and apply the following changes:

### a) Set the processing directory to the high-performance scratch

Replace:

```yaml
processing_dir: "processing"
```

with the path from `$TMPDIR/`, for example:

```yaml
processing_dir: "/tmp/3544738.1.rhel9.q/processing"
```

This tells the pipeline to copy and process input files on the high-performance local drive. Results will be copied back to the `results/` folder automatically.

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

## 7) Run the pipeline

From a conda environment with Snakemake available (e.g. `(snakemake) [mtinti@gpu-45 myRna-seq]$`):

```bash
# 11. Run the pipeline, matching the number of cores requested in qrsh
snakemake --cores 40 --configfile config_rit.yaml --use-conda
```

Sit back and watch the pipeline process the samples.

> **Tip:** Keep separate config and sample files for each experiment. This way the pipeline run is fully reproducible without needing to store intermediate files.
