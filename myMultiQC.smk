"""
MultiQC workflow with properly marked directory outputs
- Takes results directory and samples CSV file from command line
- Determines QC directories and read types from sample information
- Uses command-line options compatible with MultiQC v1.28
- Correctly marks data directories with the directory() function
- Supports both Conda and Singularity
"""
import os
import pandas as pd
import sys

# Get required parameters from command line
# Use with: snakemake --config results_dir=/path/to/results samples_csv=/path/to/samples.csv -s MultiqcOnly.smk
if "results_dir" not in config:
    raise ValueError("Please specify results_dir: snakemake --config results_dir=/path/to/results")

if "samples_csv" not in config:
    raise ValueError("Please specify samples_csv: snakemake --config samples_csv=/path/to/samples.csv")

# Validate results directory
results_dir = config['results_dir']
if not os.path.exists(results_dir):
    raise ValueError(f"Results directory not found: {results_dir}")

print(f"Using results directory: {results_dir}")

# Validate and load samples file
samples_csv = config['samples_csv']
if not os.path.exists(samples_csv):
    raise ValueError(f"Samples file not found: {samples_csv}")

print(f"Loading samples from: {samples_csv}")

try:
    samples_df = pd.read_csv(samples_csv)
    required_columns = ['sample_name', 'read_type']
    missing_columns = [col for col in required_columns if col not in samples_df.columns]
    
    if missing_columns:
        raise ValueError(f"Missing required columns in samples file: {', '.join(missing_columns)}")
        
    print(f"Loaded {len(samples_df)} samples from CSV file")
except Exception as e:
    raise ValueError(f"Error loading samples file: {str(e)}")

# Helper function to get output path
def get_output_path(*args):
    """Get a path within the MultiQC output directory"""
    return os.path.join(config['results_dir'], "MultiQC", *args)

# Find QC directories based on sample file
def find_qc_dirs_by_type():
    """Find QC directories based on sample file and categorize by read type"""
    paired_qc_dirs = []
    single_qc_dirs = []
    
    for _, row in samples_df.iterrows():
        sample_name = row['sample_name']
        read_type = row['read_type'].lower()
        
        # Check if QC directory exists for this sample
        qc_dir = os.path.join(results_dir, sample_name, "qc")
        if not os.path.isdir(qc_dir):
            print(f"Warning: QC directory not found for sample {sample_name}")
            continue
            
        # Categorize by read type from sample file
        if read_type == 'paired':
            paired_qc_dirs.append(qc_dir)
            print(f"Found paired-end sample: {sample_name}")
        elif read_type == 'single':
            single_qc_dirs.append(qc_dir)
            print(f"Found single-end sample: {sample_name}")
        else:
            print(f"Warning: Unknown read type '{read_type}' for sample {sample_name}")
    
    print(f"Found {len(paired_qc_dirs)} paired-end samples and {len(single_qc_dirs)} single-end samples")
    return paired_qc_dirs, single_qc_dirs

# Get QC directories by type
PAIRED_QC_DIRS, SINGLE_QC_DIRS = find_qc_dirs_by_type()

# Determine which reports to generate
GENERATE_PAIRED = len(PAIRED_QC_DIRS) > 0
GENERATE_SINGLE = len(SINGLE_QC_DIRS) > 0
MIXED_READ_TYPES = GENERATE_PAIRED and GENERATE_SINGLE

# Define output files based on detected read types
if MIXED_READ_TYPES:
    # If both types exist, create separate reports
    EXPECTED_OUTPUTS = [
        get_output_path("multiqc_paired_report.html"),
        get_output_path("multiqc_single_report.html")
    ]
    DATA_DIR_OUTPUTS = [
        directory(get_output_path("multiqc_paired_report_data")),
        directory(get_output_path("multiqc_single_report_data"))
    ]
elif GENERATE_PAIRED:
    # If only paired-end exists, create one report
    EXPECTED_OUTPUTS = [get_output_path("multiqc_report.html")]
    DATA_DIR_OUTPUTS = [directory(get_output_path("multiqc_report_data"))]
elif GENERATE_SINGLE:
    # If only single-end exists, create one report
    EXPECTED_OUTPUTS = [get_output_path("multiqc_report.html")]
    DATA_DIR_OUTPUTS = [directory(get_output_path("multiqc_report_data"))]
else:
    # No samples found
    EXPECTED_OUTPUTS = []
    DATA_DIR_OUTPUTS = []

# Default rule to execute everything
rule all:
    input:
        EXPECTED_OUTPUTS + DATA_DIR_OUTPUTS

# Single rule that creates the directory and generates the reports
rule generate_multiqc_reports:
    output:
        reports = EXPECTED_OUTPUTS,
        data_dirs = DATA_DIR_OUTPUTS
    params:
        multiqc_dir = get_output_path(""),
        paired_qc_dirs = " ".join(PAIRED_QC_DIRS) if PAIRED_QC_DIRS else "",
        single_qc_dirs = " ".join(SINGLE_QC_DIRS) if SINGLE_QC_DIRS else "",
        mixed_read_types = MIXED_READ_TYPES,
        generate_paired = GENERATE_PAIRED,
        generate_single = GENERATE_SINGLE,
        # Common MultiQC arguments - using version 1.28 compatible options
        common_args = (
            f"--ignore '*_unique.txt.summary' --config multiqc_config.yaml "
            f"--force -f "
        )
    conda:
        "envs/multiqc.yaml"
    singularity:
        config.get("singularity_image", "docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0")
    run:
        # Create the MultiQC directory first
        shell(f"mkdir -p {params.multiqc_dir}")
        
        if params.mixed_read_types:
            # Generate paired-end report
            paired_output = [o for o in output.reports if "paired" in o][0]
            shell(
                f"multiqc {params.paired_qc_dirs} "
                f"{params.common_args} "
                f"-o {params.multiqc_dir} "
                f"-n $(basename {paired_output}) "
                f"--title 'Paired-End Samples QC Report'"
            )

            # Generate single-end report
            single_output = [o for o in output.reports if "single" in o][0]
            shell(
                f"multiqc {params.single_qc_dirs} "
                f"{params.common_args} "
                f"-o {params.multiqc_dir} "
                f"-n $(basename {single_output}) "
                f"--title 'Single-End Samples QC Report'"
            )
        elif params.generate_paired:
            # Generate only paired-end report
            shell(
                f"multiqc {params.paired_qc_dirs} "
                f"{params.common_args} "
                f"-o {params.multiqc_dir} "
                f"-n $(basename {output.reports[0]}) "
                f"--title 'Paired-End Samples QC Report'"
            )
        elif params.generate_single:
            # Generate only single-end report
            shell(
                f"multiqc {params.single_qc_dirs} "
                f"{params.common_args} "
                f"-o {params.multiqc_dir} "
                f"-n $(basename {output.reports[0]}) "
                f"--title 'Single-End Samples QC Report'"
            )
        else:
            print("No QC directories found. No reports generated.")
            # Create empty reports to satisfy Snakemake
            for report in output.reports:
                shell(f"touch {report}")
            for data_dir in output.data_dirs:
                shell(f"mkdir -p {data_dir}")