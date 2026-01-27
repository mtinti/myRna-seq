"""
RNA-seq Snakemake Pipeline - Modular Implementation with Run Tag Support
Main workflow file that selectively includes rule modules based on sample types
"""
import os
import pandas as pd
import sys
from snakemake.workflow import workflow

# Import utility functions from lib directory
sys.path.insert(0, os.path.join(os.path.dirname(workflow.snakefile), "rules"))
from sample_utils import (
    ensure_directories_exist,
    copy_reference_files,
    load_samples,
    validate_reference_inputs,
    is_sample_completed
)

# Import configuration
configfile: "config.yaml"

# Set defaults for critical config values if not set
if 'processing_dir' not in config or config['processing_dir'] is None:
    config['processing_dir'] = "processing"
    print(f"Warning: processing_dir not set in config, using default: {config['processing_dir']}")

if 'results_dir' not in config or config['results_dir'] is None:
    config['results_dir'] = "results"
    print(f"Warning: results_dir not set in config, using default: {config['results_dir']}")

if 'reference_fasta' not in config or config['reference_fasta'] in (None, ""):
    genome_index = config.get('genome_index')
    if genome_index:
        genome_index = os.path.expanduser(genome_index)
        if genome_index.lower().endswith('.fa'):
            derived_fasta = genome_index
        else:
            derived_fasta = f"{genome_index}.fa"

        config['reference_fasta'] = derived_fasta
        print(
            "reference_fasta not provided; using FASTA derived from genome_index: "
            f"{config['reference_fasta']}"
        )

    if 'reference_fasta' not in config or config['reference_fasta'] in (None, ""):
        raise ValueError(
            "reference_fasta must be set in config.yaml or provided via --config. "
            "You may also supply genome_index pointing to the FASTA prefix (without the .fa extension)."
        )

if not str(config['reference_fasta']).lower().endswith('.fa'):
    raise ValueError(
        "reference_fasta must point to a file ending in .fa. "
        f"Received: {config['reference_fasta']}"
    )

if 'gtf_file' not in config or config['gtf_file'] is None:
    raise ValueError("gtf_file must be set in config.yaml")

if 'samples_csv' not in config or config['samples_csv'] is None:
    raise ValueError("samples_csv must be set in config.yaml")

if 'selected_samples' not in config:
    config['selected_samples'] = []

# Set global resource limits
if "max_resources" in config:
    for resource, limit in config["max_resources"].items():
        workflow.global_resources[resource] = limit

# Expand environment variables in config
for key, value in config.items():
    if isinstance(value, str) and value.startswith('$'):
        env_var = value[1:]
        if env_var in os.environ:
            config[key] = os.environ[env_var]
            print(f"Expanded environment variable in config['{key}']: {value} -> {config[key]}")
        else:
            print(f"Warning: Environment variable '{env_var}' not found. Using '{value}' as literal value.")

# Validate reference inputs
try:
    validate_reference_inputs(config)
except (FileNotFoundError, ValueError) as e:
    print(e)
    sys.exit(1)  # Exit with an error code

# Helper function to get results path (needed for is_sample_completed before including common.smk)
def get_results_path(*args):
    """Get a path within the results directory specified in config"""
    return os.path.join(config['results_dir'], *args)

def get_processing_path(*args):
    """Get a path within the processing directory specified in config"""
    return os.path.join(config['processing_dir'], *args)

# Initialize directory structure before loading samples
try:
    ensure_directories_exist(config)
    copy_reference_files(config)
except Exception as e:
    print(f"Error during initialization: {str(e)}")
    # Set default paths if there was an error
    if "processing_genome_index" not in config:
        fasta_name = os.path.splitext(os.path.basename(config["reference_fasta"]))[0]
        config["processing_genome_index"] = os.path.join(config['processing_dir'], 'reference', fasta_name)
    if "processing_annotation_source" not in config:
        config["processing_annotation_source"] = os.path.join(
            config['processing_dir'], 'reference', os.path.basename(config["gtf_file"])
        )
    if "processing_gtf_file" not in config:
        gtf_basename = os.path.basename(config["gtf_file"])
        annotation_stem, annotation_ext = os.path.splitext(gtf_basename)
        if annotation_stem.endswith(".gff"):
            annotation_stem = os.path.splitext(annotation_stem)[0]
        if annotation_ext.lower() in (".gff", ".gff3"):
            derived_gtf = f"{annotation_stem}.gtf"
        else:
            derived_gtf = gtf_basename
        config["processing_gtf_file"] = os.path.join(
            config['processing_dir'], 'reference', derived_gtf
        )
    if "processing_reference_fasta" not in config:
        config["processing_reference_fasta"] = os.path.join(
            config['processing_dir'], 'reference', os.path.basename(config["reference_fasta"]) or "reference.fa"
        )

# Load sample data - this must happen before rule parsing
SAMPLES_DF, RUN_TAG_SAMPLES, RUN_TAGS, EFFECTIVE_SAMPLES = load_samples(config)

# Determine which samples to process
all_samples = list(SAMPLES_DF.index) if not SAMPLES_DF.empty else []
completed_samples = [s for s in all_samples if is_sample_completed(config, s, SAMPLES_DF)]

# Filter out completed samples 
SAMPLES = [s for s in all_samples if s not in completed_samples]

# Create RUN_TAGS_TO_PROCESS - only run tags that have at least one sample to process
RUN_TAGS_TO_PROCESS = []
for run_tag in RUN_TAGS:
    samples = RUN_TAG_SAMPLES[run_tag]
    if any(s in SAMPLES for s in samples):
        RUN_TAGS_TO_PROCESS.append(run_tag)

# EFFECTIVE_SAMPLES_TO_PROCESS is the list of run tags and standalone samples to process
EFFECTIVE_SAMPLES_TO_PROCESS = []
for sample in SAMPLES:
    if 'run_tag' in SAMPLES_DF.columns:
        run_tag = SAMPLES_DF.loc[sample, 'run_tag']
        if pd.isna(run_tag) or run_tag == '':
            EFFECTIVE_SAMPLES_TO_PROCESS.append(sample)
        elif run_tag not in EFFECTIVE_SAMPLES_TO_PROCESS:
            EFFECTIVE_SAMPLES_TO_PROCESS.append(run_tag)
    else:
        EFFECTIVE_SAMPLES_TO_PROCESS.append(sample)

# Categorize samples by type
PAIRED_LOCAL_SAMPLES = [s for s in SAMPLES if (SAMPLES_DF.loc[s, 'read_type'] == 'paired') and (SAMPLES_DF.loc[s, 'source_type'] == 'local')]
PAIRED_FTP_SAMPLES = [s for s in SAMPLES if (SAMPLES_DF.loc[s, 'read_type'] == 'paired') and (SAMPLES_DF.loc[s, 'source_type'] == 'ftp')]
SRA_PAIRED_SAMPLES = [s for s in SAMPLES if (SAMPLES_DF.loc[s, 'read_type'] == 'paired') and (SAMPLES_DF.loc[s, 'source_type'] == 'sra_paired')]
SINGLE_LOCAL_SAMPLES = [s for s in SAMPLES if (SAMPLES_DF.loc[s, 'read_type'] == 'single') and (SAMPLES_DF.loc[s, 'source_type'] == 'local')]
SINGLE_FTP_SAMPLES = [s for s in SAMPLES if (SAMPLES_DF.loc[s, 'read_type'] == 'single') and (SAMPLES_DF.loc[s, 'source_type'] == 'ftp')]
NANOPORE_LOCAL_SAMPLES = [s for s in SAMPLES if (SAMPLES_DF.loc[s, 'read_type'] == 'nanopore') and (SAMPLES_DF.loc[s, 'source_type'] == 'local')]
NANOPORE_FTP_SAMPLES = [s for s in SAMPLES if (SAMPLES_DF.loc[s, 'read_type'] == 'nanopore') and (SAMPLES_DF.loc[s, 'source_type'] == 'ftp')]

# Print sample categorization
print(f"Paired-end, local samples: {len(PAIRED_LOCAL_SAMPLES)}")
print(f"Paired-end, FTP samples: {len(PAIRED_FTP_SAMPLES)}")
print(f"Paired-end, SRA samples: {len(SRA_PAIRED_SAMPLES)}")
print(f"Single-end, local samples: {len(SINGLE_LOCAL_SAMPLES)}")
print(f"Single-end, FTP samples: {len(SINGLE_FTP_SAMPLES)}")
print(f"Nanopore, local samples: {len(NANOPORE_LOCAL_SAMPLES)}")
print(f"Nanopore, FTP samples: {len(NANOPORE_FTP_SAMPLES)}")
print(f"Run tags to process: {len(RUN_TAGS_TO_PROCESS)}")
print(f"Effective samples to process: {len(EFFECTIVE_SAMPLES_TO_PROCESS)}")

# Cache sample availability flags for conditional module includes
HAS_LOCAL_PAIRED = bool(PAIRED_LOCAL_SAMPLES)
HAS_FTP_PAIRED = bool(PAIRED_FTP_SAMPLES)
HAS_SRA_PAIRED = bool(SRA_PAIRED_SAMPLES)
HAS_SINGLE_LOCAL = bool(SINGLE_LOCAL_SAMPLES)
HAS_SINGLE_FTP = bool(SINGLE_FTP_SAMPLES)
HAS_NANOPORE_LOCAL = bool(NANOPORE_LOCAL_SAMPLES)
HAS_NANOPORE_FTP = bool(NANOPORE_FTP_SAMPLES)
HAS_SINGLE_SHORT = HAS_SINGLE_LOCAL or HAS_SINGLE_FTP
HAS_NANOPORE = HAS_NANOPORE_LOCAL or HAS_NANOPORE_FTP
HAS_REMOTE_PAIRED = HAS_FTP_PAIRED or HAS_SRA_PAIRED

# Create a logger function that only prints at startup but not during DAG building
verbose_logging_done = False

# The onstart handler will run once before execution, but after DAG building
onstart:
    global verbose_logging_done
    
    print("RNA-seq pipeline starting...")
    
    # Print information about samples
    if completed_samples:
        print(f"Skipping {len(completed_samples)} already completed samples: {', '.join(completed_samples)}")
    
    if SAMPLES:
        print(f"Processing {len(SAMPLES)} samples: {', '.join(SAMPLES)}")
    else:
        print("Warning: No samples to process!")
    
    # Print information about run tags
    if RUN_TAGS_TO_PROCESS:
        print(f"Processing {len(RUN_TAGS_TO_PROCESS)} run tags:")
        for run_tag in RUN_TAGS_TO_PROCESS:
            print(f"  {run_tag}: {len(RUN_TAG_SAMPLES[run_tag])} samples - {', '.join(RUN_TAG_SAMPLES[run_tag])}")
    
    verbose_logging_done = True

# Include common functions
include: "rules/common.smk"

# Prepare annotation for downstream steps (convert GFF->GTF, validate feature/attribute)
include: "rules/processing/prepare_annotation.smk"

# Build Bowtie2 indexes from the staged reference FASTA
include: "rules/processing/build_bowtie_index.smk"

# Determine which run tags are non-nanopore for certain steps
NON_NANOPORE_EFFECTIVE_SAMPLES = [rt for rt in EFFECTIVE_SAMPLES_TO_PROCESS if not is_nanopore(rt)]
HAS_NON_NANOPORE_EFFECTIVE = bool(NON_NANOPORE_EFFECTIVE_SAMPLES)
HAS_PAIRED_EFFECTIVE = any(is_paired_end(run_tag) for run_tag in EFFECTIVE_SAMPLES_TO_PROCESS)
HAS_SINGLE_NON_NANO_EFFECTIVE = any(
    is_single_end(run_tag) and not is_nanopore(run_tag) for run_tag in EFFECTIVE_SAMPLES_TO_PROCESS
)
HAS_SINGLE_OR_NANO_EFFECTIVE = any(
    is_single_end(run_tag) or is_nanopore(run_tag) for run_tag in EFFECTIVE_SAMPLES_TO_PROCESS
)
HAS_PAIRED_NON_NANO_EFFECTIVE = any(
    is_paired_end(run_tag) and not is_nanopore(run_tag) for run_tag in EFFECTIVE_SAMPLES_TO_PROCESS
)

# Conditionally include rules based on sample types using cached predicates
conditional_modules = [
    ("rules/acquisition/local_paired.smk", HAS_LOCAL_PAIRED),
    ("rules/acquisition/ftp_paired.smk", HAS_FTP_PAIRED),
    ("rules/acquisition/sra_paired.smk", HAS_SRA_PAIRED),
    (
        "rules/acquisition/local_single.smk",
        HAS_SINGLE_LOCAL or HAS_NANOPORE_LOCAL,
    ),
    (
        "rules/acquisition/ftp_single.smk",
        HAS_SINGLE_FTP or HAS_NANOPORE_FTP,
    ),
    ("rules/checksum/paired_checksum.smk", HAS_LOCAL_PAIRED or HAS_FTP_PAIRED),
    ("rules/checksum/single_checksum.smk", HAS_SINGLE_SHORT or HAS_NANOPORE),
    (
        "rules/fastp/paired_fastp.smk",
        HAS_LOCAL_PAIRED or HAS_REMOTE_PAIRED,
    ),
    ("rules/fastp/single_fastp.smk", HAS_SINGLE_SHORT),
    (
        "rules/alignment/paired_align.smk",
        HAS_LOCAL_PAIRED or HAS_REMOTE_PAIRED,
    ),
    ("rules/alignment/single_align.smk", HAS_SINGLE_SHORT),
    ("rules/alignment/nanopore_align.smk", HAS_NANOPORE),
    ("rules/processing/mark_duplicates.smk", HAS_NON_NANOPORE_EFFECTIVE),
    ("rules/qc/paired_rnaseq.smk", HAS_PAIRED_NON_NANO_EFFECTIVE),
    ("rules/qc/single_rnaseq.smk", HAS_SINGLE_NON_NANO_EFFECTIVE),
    ("rules/coverage/paired_coverage.smk", HAS_PAIRED_EFFECTIVE),
    ("rules/coverage/single_coverage.smk", HAS_SINGLE_OR_NANO_EFFECTIVE),
    ("rules/counting/paired_counts.smk", HAS_PAIRED_EFFECTIVE),
    ("rules/counting/single_counts.smk", HAS_SINGLE_OR_NANO_EFFECTIVE),
    ("rules/counting/rpkm.smk", HAS_PAIRED_EFFECTIVE or HAS_SINGLE_OR_NANO_EFFECTIVE),
]

for module, should_include in conditional_modules:
    if should_include:
        include: module

# Modules that are always required
include: "rules/processing/merge_bams.smk"
include: "rules/qc/bam_qc.smk"
include: "rules/benchmarks.smk"
include: "rules/results.smk"

# Define target files outside the rule
FINAL_RESULTS = get_results_path("copy_complete_all.txt") if EFFECTIVE_SAMPLES_TO_PROCESS else []
SAMPLE_CLEANUPS = expand(get_results_path("{run_tag}/sample_cleanup_complete.txt"), run_tag=EFFECTIVE_SAMPLES_TO_PROCESS) if EFFECTIVE_SAMPLES_TO_PROCESS else []

rule all:
    input:
        # Only need to wait for results copying completion
        # Cleanup happens as part of the copy_results rule now
        get_results_path("copy_complete_all.txt") if EFFECTIVE_SAMPLES_TO_PROCESS else []
        
# Define the complete workflow steps for logging
workflow_steps = [
    "01 - Input Handling",
    "02 - Checksum Verification",
    "03 - Quality Filtering",
    "04 - Genome Alignment",
    "05 - BAM Merging by Run Tag",
    "06 - Mark Duplicates",
    "07 - BAM Quality Control",
    "08 - Coverage Tracks",
    "09 - Feature Counting",
    "10 - RPKM Calculation",
    "11 - Benchmark Analysis",
    "12 - Copy Results"
]

onsuccess:
    print("RNA-seq pipeline completed successfully!")
    print(f"Completed steps: {', '.join(workflow_steps)}")
    print(f"Results are available in: {config['results_dir']}")
    
    # Check for cleanup configuration
    if config.get("cleanup_processing", False):
        print("Processing directories have been cleaned up")
    else:
        print("Processing directories were preserved (cleanup_processing=False)")
        print(f"Intermediate files are available in: {config['processing_dir']}")

onerror:
    print("RNA-seq pipeline encountered an error.")
    print("Check the log files in the processing directory for details.")
    print("You can resume the pipeline with --keep-going to continue from where it left off.")
