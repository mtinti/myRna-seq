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
    validate_genome_index,
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

if 'genome_index' not in config or config['genome_index'] is None:
    raise ValueError("genome_index must be set in config.yaml")

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

# Validate reference genome index
try:
    validate_genome_index(config)
except FileNotFoundError as e:
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
        config["processing_genome_index"] = os.path.join(config['processing_dir'], 'reference', os.path.basename(config["genome_index"]))
    if "processing_gtf_file" not in config:
        config["processing_gtf_file"] = os.path.join(config['processing_dir'], 'reference', os.path.basename(config["gtf_file"]))

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
SINGLE_LOCAL_SAMPLES = [s for s in SAMPLES if (SAMPLES_DF.loc[s, 'read_type'] == 'single') and (SAMPLES_DF.loc[s, 'source_type'] == 'local')]
SINGLE_FTP_SAMPLES = [s for s in SAMPLES if (SAMPLES_DF.loc[s, 'read_type'] == 'single') and (SAMPLES_DF.loc[s, 'source_type'] == 'ftp')]

# Print sample categorization
print(f"Paired-end, local samples: {len(PAIRED_LOCAL_SAMPLES)}")
print(f"Paired-end, FTP samples: {len(PAIRED_FTP_SAMPLES)}")
print(f"Single-end, local samples: {len(SINGLE_LOCAL_SAMPLES)}")
print(f"Single-end, FTP samples: {len(SINGLE_FTP_SAMPLES)}")
print(f"Run tags to process: {len(RUN_TAGS_TO_PROCESS)}")
print(f"Effective samples to process: {len(EFFECTIVE_SAMPLES_TO_PROCESS)}")

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

# Conditionally include rules based on sample types
# Acquisition rules
if PAIRED_LOCAL_SAMPLES:
    include: "rules/acquisition/local_paired.smk"

if PAIRED_FTP_SAMPLES:
    include: "rules/acquisition/ftp_paired.smk"

if SINGLE_LOCAL_SAMPLES:
    include: "rules/acquisition/local_single.smk"

if SINGLE_FTP_SAMPLES:
    include: "rules/acquisition/ftp_single.smk"

# Checksum verification rules
if len(PAIRED_LOCAL_SAMPLES) + len(PAIRED_FTP_SAMPLES) > 0:
    include: "rules/checksum/paired_checksum.smk"

if len(SINGLE_LOCAL_SAMPLES) + len(SINGLE_FTP_SAMPLES) > 0:
    include: "rules/checksum/single_checksum.smk"

# Quality filtering rules
if len(PAIRED_LOCAL_SAMPLES) + len(PAIRED_FTP_SAMPLES) > 0:
    include: "rules/fastp/paired_fastp.smk"

if len(SINGLE_LOCAL_SAMPLES) + len(SINGLE_FTP_SAMPLES) > 0:
    include: "rules/fastp/single_fastp.smk"

# Alignment rules
if len(PAIRED_LOCAL_SAMPLES) + len(PAIRED_FTP_SAMPLES) > 0:
    include: "rules/alignment/paired_align.smk"

if len(SINGLE_LOCAL_SAMPLES) + len(SINGLE_FTP_SAMPLES) > 0:
    include: "rules/alignment/single_align.smk"

# BAM merging rule (new module)
include: "rules/processing/merge_bams.smk"

# Common processing rules
include: "rules/processing/mark_duplicates.smk"

# QC rules
include: "rules/qc/bam_qc.smk"

if any(is_paired_end(run_tag) for run_tag in EFFECTIVE_SAMPLES_TO_PROCESS):
    include: "rules/qc/paired_rnaseq.smk"

if any(is_single_end(run_tag) for run_tag in EFFECTIVE_SAMPLES_TO_PROCESS):
    include: "rules/qc/single_rnaseq.smk"

# Coverage track generation rules
if any(is_paired_end(run_tag) for run_tag in EFFECTIVE_SAMPLES_TO_PROCESS):
    include: "rules/coverage/paired_coverage.smk"

if any(is_single_end(run_tag) for run_tag in EFFECTIVE_SAMPLES_TO_PROCESS):
    include: "rules/coverage/single_coverage.smk"

# Feature counting rules
if any(is_paired_end(run_tag) for run_tag in EFFECTIVE_SAMPLES_TO_PROCESS):
    include: "rules/counting/paired_counts.smk"

if any(is_single_end(run_tag) for run_tag in EFFECTIVE_SAMPLES_TO_PROCESS):
    include: "rules/counting/single_counts.smk"

# Benchmarking rules
include: "rules/benchmarks.smk"

# Results handling rules
include: "rules/results.smk"

# Define final target rule
rule all:
    input:
        # First process all original samples through alignment
        acquisition_flags = [get_acquisition_flag(sample) for sample in SAMPLES],
        checksum_flags = [get_checksum_flag(sample) for sample in SAMPLES],
        fastp_flags = [get_fastp_flag(sample) for sample in SAMPLES],
        alignment_flags = [get_alignment_flag(sample) for sample in SAMPLES],
        
        # Then merge and process all effective samples (run tags or standalone samples)
        merge_flags = [get_merge_flag(run_tag) for run_tag in EFFECTIVE_SAMPLES_TO_PROCESS],
        markduplicates_flags = [get_markduplicates_flag(run_tag) for run_tag in EFFECTIVE_SAMPLES_TO_PROCESS],
        qc_flags = [get_qc_flag(run_tag) for run_tag in EFFECTIVE_SAMPLES_TO_PROCESS],
        coverage_flags = [get_coverage_flag(run_tag) for run_tag in EFFECTIVE_SAMPLES_TO_PROCESS],
        featurecounts_flags = [get_featurecounts_flag(run_tag) for run_tag in EFFECTIVE_SAMPLES_TO_PROCESS],
        
        # Final results flags
        results_complete = get_results_path("copy_complete_all.txt") if EFFECTIVE_SAMPLES_TO_PROCESS else [],
        cleanup_complete = get_results_path("cleanup_complete.txt") if config.get("cleanup_processing", False) and EFFECTIVE_SAMPLES_TO_PROCESS else []
        
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
    "10 - Benchmark Analysis",
    "11 - Copy Results"
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