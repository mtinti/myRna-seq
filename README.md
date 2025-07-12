# RNA-seq Snakemake Pipeline

A comprehensive RNA-seq analysis pipeline built with Snakemake that supports multiple data types including paired-end, single-end, interleaved, and nanopore sequencing data.

## Features

- **Multi-format support**: Paired-end, single-end, interleaved, and nanopore data
- **Flexible input sources**: Local files and FTP downloads
- **Quality control**: Comprehensive QC with fastp, samtools, and Qualimap
- **Run tag grouping**: Combine multiple samples for joint analysis
- **Checksum verification**: Ensure data integrity
- **Containerized**: Singularity support
- **Benchmarking**: Detailed performance metrics
- **Modular design**: Rules are automatically included based on your data
- **On-demand directories**: Sample folders are created only when processing starts and cleaned up when finished

## Pipeline Architecture

```
Input Data (4 types)
    ↓
Acquisition & Checksum Verification
    ↓
Quality Filtering* (fastp)
    ↓
Alignment (bowtie2/minimap2)
    ↓
BAM Merging by Run Tag
    ↓
Mark Duplicates* (Picard)
    ↓
Quality Control (samtools + Qualimap*)
    ↓
Coverage Tracks (bamCoverage)
    ↓
Feature Counting (featureCounts)
    ↓
Results & Cleanup

* Step skipped for nanopore data
```

### Data Type Processing:
- **Paired/Single/Interleaved**: Full pipeline with bowtie2
- **Nanopore**: Direct minimap2 alignment, skips QC/duplicate steps

## Supported Data Types

### 1. Paired-end Illumina Data
- **read_type**: `paired`
- **Description**: Standard paired-end sequencing with separate R1 and R2 files
- **Required fields**: `file_path_1` (R1), `file_path_2` (R2)
- **Processing**: Quality filtering with fastp → Alignment with bowtie2 → Standard QC

### 2. Single-end Illumina Data
- **read_type**: `single`
- **Description**: Single-end sequencing with one file per sample
- **Required fields**: `file_path_1` only
- **Processing**: Quality filtering with fastp → Alignment with bowtie2 → Standard QC

### 3. Interleaved Paired-end Data
- **read_type**: `interleaved`
- **Description**: Paired-end data stored in a single file with alternating R1/R2 reads
- **Required fields**: `file_path_1` (interleaved file)
- **Processing**: File splitting → Quality filtering with fastp → Alignment with bowtie2 → Standard QC
- **Note**: The pipeline automatically splits interleaved files into separate R1/R2 files before processing

### 4. Nanopore Long-read Data
- **read_type**: `nanopore`
- **Description**: Oxford Nanopore long-read sequencing data
- **Required fields**: `file_path_1` only
- **Processing**: Direct alignment with minimap2 → Coverage analysis (skips quality filtering and duplicate marking)
- **Special considerations**:
  - Uses minimap2 with splice-aware alignment (`-ax splice`)
  - Skips fastp quality filtering (nanopore reads handled differently)
  - Skips Picard MarkDuplicates (not applicable to long reads)
  - Skips Qualimap RNA-seq QC (optimized for short reads)

## Quick Start

### 1. Installation

```bash
# Clone the repository
git clone <repository-url>
cd rna-seq-pipeline

# Install Snakemake (if not already installed)
conda install -c conda-forge -c bioconda snakemake
```

### 2. Prepare your sample sheet

Create a `samples.csv` file with your sample information:

```csv
sample_name,read_type,source_type,file_path_1,file_path_2,checksum_1,checksum_2,run_tag
sample1,paired,local,/data/sample1_R1.fastq.gz,/data/sample1_R2.fastq.gz,abc123,def456,experiment_A
sample2,single,ftp,ftp://server.com/sample2.fastq.gz,,xyz789,,experiment_B
sample3,interleaved,local,/data/sample3_interleaved.fastq.gz,,mno345,,experiment_A
sample4,nanopore,local,/data/sample4_nanopore.fastq.gz,,pqr678,,long_read_set
```

### 3. Configure the pipeline

Edit `config.yaml` to specify your reference genome and processing parameters:

```yaml
# Core directories
processing_dir: "processing"
results_dir: "results"

# Sample information
samples_csv: "samples.csv"

# Reference genome
genome_index: "reference/genome/genome_index"
gtf_file: "reference/genome/annotation.gtf"

# Processing parameters
cores_align: 8
feature_type: "CDS"
```

### 4. Run the pipeline

```bash
# Dry run to check the workflow
snakemake -n

# Run the pipeline
snakemake --cores 8

# Run with Singularity container
snakemake --cores 8 --use-singularity
```

## Testing the Pipeline

You can test the pipeline with the included test dataset:

```bash
snakemake --cores 10 --use-singularity --config \
processing_dir="tests/test_counts/new_branch/processing" \
results_dir="tests/test_counts/new_branch/results" \
benchmark_dir="tests/test_counts/new_branch/benchmarks" \
genome_index="tests/test_counts/genome/random_genome" \
gtf_file="tests/test_counts/genome/annotation.gtf" \
samples_csv="test_samples_counts.csv"
```

Expected output from the test dataset:
```
Features
---------------------------
gene1: 10 reads total
gene2: 8 reads total
gene3: 2 reads total
---------------------------
```

## Sample CSV Format

Your `samples.csv` file should include these columns:

| Column | Required | Description |
|--------|----------|-------------|
| `sample_name` | Yes | Unique identifier for each sample |
| `read_type` | Yes | One of: `paired`, `single`, `interleaved`, `nanopore` |
| `source_type` | Yes | `local`, `ftp`, or `sra_paired` |
| `file_path_1` | Yes | Path to first/only FASTQ file |
| `file_path_2` | Conditional | Path to second FASTQ file (required for `paired` only) |
| `checksum_1` | Optional | MD5 checksum for file_path_1 |
| `checksum_2` | Optional | MD5 checksum for file_path_2 |
| `run_tag` | Optional | Group samples for combined analysis |

### Example CSV entries:

```csv
sample_name,read_type,source_type,file_path_1,file_path_2,checksum_1,checksum_2,run_tag
illumina_pe_01,paired,local,/data/sample1_R1.fastq.gz,/data/sample1_R2.fastq.gz,abc123,def456,experiment_A
illumina_se_01,single,ftp,ftp://server.com/sample2.fastq.gz,,xyz789,,experiment_B
interleaved_01,interleaved,local,/data/sample3_interleaved.fastq.gz,,mno345,,experiment_A
nanopore_01,nanopore,local,/data/sample4_nanopore.fastq.gz,,pqr678,,long_read_set
```

## Using SRA Data

The pipeline can download reads directly from NCBI SRA. Set `source_type` to `sra_paired` and specify the accession in `file_path_1`:

```csv
sample_name,read_type,source_type,file_path_1,file_path_2,checksum_1,checksum_2
MySRA,paired,sra_paired,SRR12345678,,,
```

SRA reads are downloaded with `fastq-dump`, reformatted, and gzipped before entering the standard workflow.

## Processing Differences by Data Type

### Quality Control Steps

| Step | Paired/Single | Interleaved | Nanopore |
|------|---------------|-------------|----------|
| File splitting | No | ✅ Yes | No |
| fastp filtering | ✅ Yes | ✅ Yes | ❌ No |
| Alignment tool | bowtie2 | bowtie2 | minimap2 |
| Mark duplicates | ✅ Yes | ✅ Yes | ❌ No |
| Qualimap BAM QC | ✅ Yes | ✅ Yes | ❌ No |
| Qualimap RNA-seq QC | ✅ Yes | ✅ Yes | ❌ No |

### Alignment Parameters

**Illumina data (paired/single/interleaved):**
- Tool: bowtie2
- Options: `--very-sensitive-local`
- Read groups: Automatically added for Picard compatibility

**Nanopore data:**
- Tool: minimap2
- Preset: `map-ont` (configurable via `nanopore_preset`)
- Options: `-ax splice -uf -k14 -G 10000` (splice-aware alignment)
- Read groups: Platform set to "ONT"

## Configuration

### Core Configuration (`config.yaml`)

```yaml
# Core directories
processing_dir: "processing"
results_dir: "results"

# Sample information
samples_csv: "samples.csv"
selected_samples: []  # Empty list means all samples are processed

# Reference genome
genome_index: "reference/genome/genome_index"
gtf_file: "reference/genome/annotation.gtf"

# Processing parameters
cores_align: 8
cores_coverage: 8
cores_featurecounts: 8
cores_fastp: 8
cores_default: 1
feature_type: "CDS"

# Nanopore-specific settings
nanopore_preset: "map-ont"        # Minimap2 preset
cores_nanopore_align: 8           # CPU cores for nanopore alignment
mem_nanopore_align: 16000         # Memory for nanopore alignment (MB)

# Mark Duplicates parameters
remove_duplicates: false  # Set to true to remove duplicates, false to only mark

# QC parameters
qualimap_memory: "10G"

# Coverage track parameters
coverage_bin_size: 10
coverage_normalize: "RPKM"
min_mapping_quality: 2

# File handling preferences
copy_fastq: false        # Whether to copy FASTQ files to results
copy_bam: false          # Whether to copy BAM files to results
copy_benchmarks: true    # Whether to copy benchmark files
copy_logs: true          # Whether to copy log files
cleanup_processing: false  # Remove processing directory and reference after completion

# Container support
singularity_image: "/path/to/container.sif"
```

### Resource Configuration

```yaml
# Resource limits
max_resources:
  network: 3    # Maximum concurrent FTP downloads
  io: 1000      # Maximum concurrent IO operations
```

## Output Structure

```
results/
├── {sample_or_run_tag}/
│   ├── {sample}_all.bw                    # Coverage track (all reads)
│   ├── {sample}_unique.bw                 # Coverage track (unique reads)
│   ├── {sample}_counts_*_all.txt          # Feature counts (all reads)
│   ├── {sample}_counts_*_unique.txt       # Feature counts (unique reads)
│   ├── qc/                                # Quality control reports
│   │   ├── fastp/                         # fastp reports (if applicable)
│   │   ├── flagstat/                      # samtools flagstat
│   │   ├── stats/                         # samtools stats
│   │   ├── {sample}_qualimap_bam/         # Qualimap BAM QC
│   │   ├── {sample}_qualimap_rnaseq/      # Qualimap RNA-seq QC (if applicable)
│   │   ├── markduplicates/                # Picard metrics (if applicable)
│   │   └── feature_counts/                # featureCounts summaries
│   ├── benchmarks/                        # Performance benchmarks (optional)
│   └── logs/                              # Log files (optional)
└── copy_complete_all.txt                  # Pipeline completion marker
```

If `cleanup_processing` is enabled, the pipeline removes each sample's processing
directory after the results are copied. The `processing/reference` directory is
also removed at the end of the run.

## Run Tags and Sample Grouping

The `run_tag` column allows you to group samples for combined analysis:

```csv
sample_name,read_type,source_type,file_path_1,run_tag
sample_A1,paired,local,/data/A1_R1.fastq.gz,experiment_1
sample_A2,paired,local,/data/A2_R1.fastq.gz,experiment_1
sample_B1,nanopore,local,/data/B1.fastq.gz,long_reads
```

**Important notes:**
- All samples within a run tag group must have the same `read_type`
- BAM files are merged by run tag before downstream analysis
- Final outputs use the run tag name, not individual sample names
- If no `run_tag` is specified, each sample is processed independently

## Advanced Usage

### Running Specific Samples

```bash
# Process only specific samples
snakemake --cores 8 --config selected_samples="['sample1','sample2']"
```

### Container Usage

```bash
# With Singularity
snakemake --cores 8 --use-singularity
```

### Resource Management

```bash
# Limit concurrent FTP downloads
snakemake --cores 8 --resources network=2
```

### Resource Management

```bash
# Limit concurrent FTP downloads
snakemake --cores 8 --resources network=2
```