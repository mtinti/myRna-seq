"""
Common functions and variables for RNA-seq pipeline
Helper functions that work with SAMPLES_DF, RUN_TAGS, and other global variables
"""
import os

# Helper function to get processing directory path
def get_processing_path(*args):
    """Get a path within the processing directory specified in config"""
    return os.path.join(config['processing_dir'], *args)

# Helper function to get results directory path
def get_results_path(*args):
    """Get a path within the results directory specified in config"""
    return os.path.join(config['results_dir'], *args)

# Helper functions for run tag information
def get_run_tag(sample):
    """Get the run tag for a sample, or the sample name if no run tag exists"""
    try:
        if 'run_tag' not in SAMPLES_DF.columns:
            return sample
        
        run_tag = SAMPLES_DF.loc[sample, 'run_tag']
        if pd.isna(run_tag) or run_tag == '':
            return sample
        return run_tag
    except Exception as e:
        print(f"Error getting run tag for sample {sample}: {str(e)}")
        return sample

def is_run_tag(name):
    """Check if a name is a run tag or an original sample name"""
    return name in RUN_TAGS

# Helper functions for sample information
def get_read_type(sample):
    """Return 'paired' or 'single' for a sample"""
    try:
        return SAMPLES_DF.loc[sample, 'read_type']
    except Exception as e:
        print(f"Error getting read type for sample {sample}: {str(e)}")
        # Default to paired
        return 'paired'

def get_source_type(sample):
    """Return 'local' or 'ftp' for a sample"""
    try:
        return SAMPLES_DF.loc[sample, 'source_type']
    except Exception as e:
        print(f"Error getting source type for sample {sample}: {str(e)}")
        # Default to local
        return 'local'

def get_run_tag_read_type(run_tag):
    """Get the read type for a run tag by checking its samples"""
    if run_tag not in RUN_TAG_SAMPLES:
        # If not a run tag or not in mapping, return the sample's own read type
        return get_read_type(run_tag)
    
    # For run tags, check the read type of the first sample
    samples = RUN_TAG_SAMPLES[run_tag]
    if not samples:
        return 'paired'  # Default if no samples
    
    return get_read_type(samples[0])

def is_paired_end(sample_or_run_tag):
    """Check if sample or run tag is paired-end"""
    if is_run_tag(sample_or_run_tag):
        return get_run_tag_read_type(sample_or_run_tag) == 'paired'
    return get_read_type(sample_or_run_tag) == 'paired'

def is_single_end(sample_or_run_tag):
    """Check if sample or run tag is single-end"""
    if is_run_tag(sample_or_run_tag):
        return get_run_tag_read_type(sample_or_run_tag) == 'single'
    return get_read_type(sample_or_run_tag) == 'single'

def is_local_source(sample):
    """Check if sample source is local"""
    return get_source_type(sample) == 'local'

def is_ftp_source(sample):
    """Check if sample source is FTP"""
    return get_source_type(sample) == 'ftp'

# Helper functions for accessing input files
def get_fastq_r1(wildcards):
    """Get the R1 fastq file for a paired-end sample"""
    sample = wildcards.sample
    return get_processing_path(f"{sample}/{sample}.1.fastq.gz")

def get_fastq_r2(wildcards):
    """Get the R2 fastq file for a paired-end sample"""
    sample = wildcards.sample
    return get_processing_path(f"{sample}/{sample}.2.fastq.gz")

def get_fastq_single(wildcards):
    """Get the fastq file for a single-end sample"""
    sample = wildcards.sample
    return get_processing_path(f"{sample}/{sample}.fastq.gz")

# Helper functions for accessing cleaned FASTQ files (after fastp)
def get_cleaned_fastq_r1(wildcards):
    """Get the cleaned R1 fastq file for a paired-end sample"""
    sample = wildcards.sample
    return get_processing_path(f"{sample}/{sample}.1.cleaned.fastq.gz")

def get_cleaned_fastq_r2(wildcards):
    """Get the cleaned R2 fastq file for a paired-end sample"""
    sample = wildcards.sample
    return get_processing_path(f"{sample}/{sample}.2.cleaned.fastq.gz")

def get_cleaned_fastq_single(wildcards):
    """Get the cleaned fastq file for a single-end sample"""
    sample = wildcards.sample
    return get_processing_path(f"{sample}/{sample}.cleaned.fastq.gz")

# Helper function for BAM file access
def get_aligned_bam(wildcards):
    """Get the aligned BAM file for a sample"""
    # Handle both direct wildcards object and dictionary input
    if isinstance(wildcards, dict):
        sample = wildcards.get("sample")
    else:
        sample = wildcards.sample
    return get_processing_path(f"{sample}/{sample}.bam")

# Helper function for accessing the Picard BAM file
def get_picard_bam(wildcards):
    """Get the BAM file after processing with Picard MarkDuplicates"""
    run_tag = wildcards.run_tag if hasattr(wildcards, 'run_tag') else wildcards.sample
    return get_processing_path(f"{run_tag}/{run_tag}.picard.bam")

# Helper functions for flag file paths with specific types
def get_acquisition_flag(sample):
    """Get the appropriate acquisition flag based on sample type"""
    if is_paired_end(sample):
        return get_processing_path(f"{sample}/acquisition_paired_complete.flag")
    else:
        return get_processing_path(f"{sample}/acquisition_single_complete.flag")

def get_checksum_flag(sample):
    """Get the appropriate checksum flag based on sample type"""
    if is_paired_end(sample):
        return get_processing_path(f"{sample}/checksums_paired_verified.flag")
    else:
        return get_processing_path(f"{sample}/checksums_single_verified.flag")

def get_fastp_flag(sample):
    """Get the appropriate fastp flag based on sample type"""
    if is_paired_end(sample):
        return get_processing_path(f"{sample}/fastp_paired_complete.flag")
    else:
        return get_processing_path(f"{sample}/fastp_single_complete.flag")

def get_alignment_flag(sample):
    """Get the appropriate alignment flag based on sample type"""
    if is_paired_end(sample):
        return get_processing_path(f"{sample}/alignment_paired_complete.flag")
    else:
        return get_processing_path(f"{sample}/alignment_single_complete.flag")

def get_merge_flag(run_tag):
    """Get the merge flag for a run tag"""
    return get_processing_path(f"merged/{run_tag}/merge_complete.flag")

def get_qc_flag(run_tag):
    """Get the appropriate QC flag based on run tag type"""
    if is_paired_end(run_tag):
        return get_processing_path(f"{run_tag}/qc_paired_complete.flag")
    else:
        return get_processing_path(f"{run_tag}/qc_single_complete.flag")

def get_coverage_flag(run_tag):
    """Get the appropriate coverage flag based on run tag type"""
    if is_paired_end(run_tag):
        return get_processing_path(f"{run_tag}/coverage_paired_complete.flag")
    else:
        return get_processing_path(f"{run_tag}/coverage_single_complete.flag")

def get_featurecounts_flag(run_tag):
    """Get the appropriate feature counting flag based on run tag type"""
    if is_paired_end(run_tag):
        return get_processing_path(f"{run_tag}/featurecounts_paired_complete.flag")
    else:
        return get_processing_path(f"{run_tag}/featurecounts_single_complete.flag")

def get_markduplicates_flag(run_tag):
    """Get the mark duplicates flag for a run tag"""
    return get_processing_path(f"{run_tag}/markduplicates_complete.flag")