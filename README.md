# 🧬 RNA-seq Processing Pipeline

A flexible and modular Snakemake workflow for processing RNA-seq data from trypanosomatids, supporting both paired-end and single-end reads from local and FTP sources.

## 📋 Overview


This pipeline was specifically designed for analysis of trypanosomatid RNA-seq data (such as *Trypanosoma brucei*, *Leishmania* spp. and related organisms). The unique genomic features of these organisms—including their intron-less genes and highly repetitive genomic regions—have driven tool selection and parameter optimization throughout the workflow.

Currently implemented in our laboratory for RNA-seq data processing prior to differential expression analysis, this pipeline handles the critical steps from raw data to quantified gene expression with parameters optimized for trypanosomatid data.

The workflow includes:

- 📥 **Acquisition**: Downloads or copies data from FTP or local sources
- ✓ **Validation**: Verifies file integrity with MD5 checksums
- 🧹 **Quality Control**: Filters low-quality reads with fastp
- 🔍 **Alignment**: Maps reads to a reference genome with Bowtie2
- 🔄 **Deduplication**: Marks or removes duplicates with Picard
- 📊 **Analysis**: Generates coverage tracks, counts reads in features
- 📈 **QC Reports**: Produces comprehensive quality metrics
- 🗃️ **Results**: Organizes outputs in a clean directory structure

## 🧪 Design Considerations for Trypanosomatid Data

The pipeline incorporates several specific optimizations for working with trypanosomatid genomic data:

1. **Intron-less Gene Handling**: Since trypanosomatids generally lack introns, the feature counting parameters are optimized for CDS-based quantification rather than spliced transcript analysis. This approach also addresses the typically poor annotation of UTRs in these organisms.

2. **Repetitive Region Management**: The pipeline provides both all-mapped and uniquely-mapped read quantification to help address the challenges of repetitive genes common in trypanosomatid genomes.

These specialized optimizations make this pipeline particularly effective for research on *T. brucei*, *T. cruzi*, *Leishmania* species, and other kinetoplastids, while still maintaining flexibility for use with other data types.

## 🔧 Pipeline Architecture

```
Input Data (Local/FTP)
      ↓
    QC FASTQ
      ↓
   Alignment
      ↓
Mark Duplicates
      ↓    
      ├─────┬─────┐
      ↓     ↓     ↓
Coverage  QC BAM  Feature Counting
      ↓     ↓     ↓
      └─────┴─────┘
      ↓
  Final Results
```



## ⚙️ Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/username/rnaseq-pipeline.git
   cd rnaseq-pipeline
   ```

2. Create a conda environment:
   ```bash
   conda env create -f environment.yaml
   conda activate rnaseq
   ```

3. Container support:
   
   The pipeline can be run with Singularity for improved reproducibility and dependency management.
   
   **Using the pre-built Docker image**
   
   A Docker image is available at Docker Hub: `mtinti/rna_seq`
   
   To convert the Docker image to a Singularity SIF file:
   
   ```bash
   # Pull from Docker Hub and create a SIF file
   singularity pull docker://mtinti/rna_seq
   
   # Or with a specific version
   singularity pull docker://mtinti/rna_seq:latest
   
   # You can also build with more options
   singularity build --sandbox rna_seq/ docker://mtinti/rna_seq
   singularity build rna_seq.sif rna_seq/
   ```
   
   The resulting SIF file can be used with the pipeline by setting `singularity_image` in the config file or at runtime.

## 🚀 Usage

### Basic Usage

Run the pipeline with default settings:

```bash
snakemake --cores 8 --use-conda
```

### Advanced Usage

Run with custom parameters:

```bash
snakemake --cores 8 --use-conda \
  --config processing_dir=/path/to/processing \
  results_dir=/path/to/results \
  samples_csv=/path/to/samples.csv
```

Run with Singularity:

```bash
# Using a local SIF file
snakemake --cores 8 --use-singularity \
  --config singularity_image="/path/to/rna_seq.sif"
```

### Command Line Arguments

The pipeline supports several command line arguments through Snakemake's `--config` parameter:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `processing_dir` | Directory for intermediate files | `processing` |
| `results_dir` | Directory for final results | `results` |
| `samples_csv` | Sample information CSV | `samples.csv` |
| `selected_samples` | Comma-separated list of specific samples to process | All samples |
| `genome_index` | Path to Bowtie2 genome index | `reference/genome` |
| `gtf_file` | Path to GTF file | `reference/genes.gtf` |

Example:

```bash
snakemake --cores 8 --use-conda \
  --config processing_dir=/scratch/user/rnaseq_work \
           results_dir=/home/user/rnaseq_results \
           samples_csv=project_samples.csv \
           selected_samples=sample1,sample2,sample5
```

## 📋 Sample File Format

The pipeline requires a CSV file with sample information. Example format:

```csv
sample_name,read_type,source_type,file_path_1,file_path_2,checksum_1,checksum_2
sample1,paired,local,R1.fastq.gz,R2.fastq.gz,5f363e2a59c7,a9c8e0d1b3f5
sample2,single,local,fastq.gz,,7d8f9e2a3b1c,
```

Required columns:
- `sample_name`: Unique sample identifier
- `read_type`: Either `paired` or `single`
- `source_type`: Either `local` or `ftp`
- `file_path_1`: Path to R1 file (or the only file for single-end)
- `file_path_2`: Path to R2 file (only for paired-end, leave empty for single-end)
- `checksum_1`: MD5 checksum for the first file (or the only file for single-end)
- `checksum_2`: MD5 checksum for the second file (only for paired-end, leave empty for single-end)

## 🔍 Output Structure

The pipeline organizes results in a clean directory structure:

```
results/
├── sample1/
│   ├── sample1_all.bw                     # Coverage track (all reads)
│   ├── sample1_unique.bw                  # Coverage track (unique reads)
│   ├── sample1_counts_paired_all.txt      # Read counts (all reads)
│   ├── sample1_counts_paired_unique.txt   # Read counts (unique reads)
│   ├── merged_benchmarks.txt              # Performance metrics
│   ├── qc/
│   │   ├── flagstat/                      # Alignment statistics
│   │   ├── stats/                         # Detailed BAM statistics
│   │   ├── sample1_qualimap_bam/          # BAM quality metrics
│   │   ├── sample1_qualimap_rnaseq/       # RNA-seq specific metrics
│   │   ├── fastp/                         # Read quality reports
│   │   ├── bowtie2/                       # Alignment reports
│   │   ├── markduplicates/                # Duplication metrics
│   │   └── feature_counts/                # Counting summaries
│   └── [optional: BAM files, logs, etc.]
├── sample2/
│   └── ...
└── copy_complete_all.txt                  # Processing completion flag
```

## 🛠 Configuration

Edit `config.yaml` to customize pipeline behavior:

```yaml
# Core settings
processing_dir: "processing"           # Intermediate files location
results_dir: "results"                 # Final output location
samples_csv: "samples.csv"             # Sample information

# Resources
cores_align: 8                         # Cores for alignment
cores_coverage: 8                      # Cores for coverage generation
cores_featurecounts: 8                 # Cores for counting
max_cores: 8                           # Maximum cores to use

# Parameters
feature_type: "CDS"                    # Feature type for counting (CDS recommended for trypanosomatids)
remove_duplicates: false               # Mark or remove duplicates
coverage_normalize: "RPKM"             # Coverage normalization method
min_mapping_quality: 2                 # Min quality for unique reads

# Trypanosomatid-specific settings
# For repetitive genes, use both 'all' and 'unique' counting approaches

# Output control
copy_bam: false                        # Copy BAM files to results
copy_benchmarks: true                  # Copy benchmark files
copy_logs: true                        # Copy log files
cleanup_processing: false              # Remove processing files after completion
```

## 🧪 Testing

Run a test with the included test data:

```bash
snakemake --cores 2 --use-conda -p --configfile test_config.yaml
```


## 🙏 Acknowledgments

This pipeline uses several excellent open-source tools:
- [Snakemake](https://snakemake.github.io)
- [fastp](https://github.com/OpenGene/fastp)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2)
- [Samtools](http://www.htslib.org)
- [Picard](https://broadinstitute.github.io/picard)
- [Qualimap](http://qualimap.conesalab.org)
- [deepTools](https://deeptools.readthedocs.io)
- [Subread/featureCounts](http://subread.sourceforge.net)
