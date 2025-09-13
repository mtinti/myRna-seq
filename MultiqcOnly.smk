"""
Simplified MultiQC workflow for RNA-seq pipeline results
- No complex file fixing or renaming
- Just parses QC folders directly
- Creates separate reports for paired-end and single-end samples
"""
import os
import pandas as pd
from pathlib import Path

# Import configuration
configfile: "config.yaml"


# Check if a custom MultiQC config was specified
if 'multiqc_config' in config:
    multiqc_config_file = config['multiqc_config']
    use_custom_multiqc_config = os.path.exists(multiqc_config_file)
    if use_custom_multiqc_config:
        print(f"Using custom MultiQC config: {multiqc_config_file}")
    else:
        print(f"Warning: Specified MultiQC config '{multiqc_config_file}' does not exist. Will generate default config.")
        use_custom_multiqc_config = False
else:
    use_custom_multiqc_config = False
    multiqc_config_file = None
    print("No custom MultiQC config specified. Will generate default config.")



# Helper function to get results directory path
def get_results_path(*args):
    """Get a path within the results directory specified in config"""
    return os.path.join(config['results_dir'], *args)


try:
    results_dir = config['results_dir']
    all_samples = [d for d in os.listdir(results_dir) 
                  if os.path.isdir(os.path.join(results_dir, d)) and 
                  os.path.exists(os.path.join(results_dir, d, "qc"))]
    print(f"Detected {len(all_samples)} samples from results directory")
except Exception as e:
    print(f"Error detecting samples from results directory: {str(e)}")
    all_samples = []

print(f"Processing {len(all_samples)} samples: {', '.join(all_samples)}")

# Helper functions for sample information
def get_read_type(sample):
    """Return 'paired' or 'single' for a sample"""
    try:
        # Detect based on QC file patterns
        sample_qc_dir = os.path.join(config['results_dir'], sample, "qc")
        # Look for bowtie2 paired/single stats file
        if os.path.exists(os.path.join(sample_qc_dir, "bowtie2", f"{sample}.bowtie2_paired_stats.txt")):
            return 'paired'
        elif os.path.exists(os.path.join(sample_qc_dir, "fastp", f"{sample}.fastp.paired.html")):
            return 'paired'
        elif os.path.exists(os.path.join(sample_qc_dir, "feature_counts", f"{sample}_counts_paired_all.txt.summary")):
            return 'paired'
        # Default to single-end if we can't find paired-end indicators
        return 'single'
    except Exception as e:
        print(f"Error determining read type for sample {sample}: {str(e)}")
        # Default to single-end to be safe
        return 'single'

def is_paired_end(sample):
    """Check if sample is paired-end"""
    return get_read_type(sample) == 'paired'

def is_single_end(sample):
    """Check if sample is single-end"""
    return get_read_type(sample) == 'single'

# Function to filter samples by read type
def get_paired_samples():
    """Get list of paired-end samples"""
    paired = [s for s in all_samples if is_paired_end(s)]
    print(f"Found {len(paired)} paired-end samples: {', '.join(paired) if paired else 'None'}")
    return paired

def get_single_samples():
    """Get list of single-end samples"""
    single = [s for s in all_samples if is_single_end(s)]
    print(f"Found {len(single)} single-end samples: {', '.join(single) if single else 'None'}")
    return single

# Evaluate sample lists once
paired_samples = get_paired_samples()
single_samples = get_single_samples()

# Default rule to execute everything
rule all:
    input:
        get_results_path("MultiQC/multiqc_reports_complete")

# Rule to create MultiQC directory
rule create_multiqc_dir:
    output:
        multiqc_dir = directory(get_results_path("MultiQC")),
        log_dir = directory(get_results_path("MultiQC/logs")),
        init_log = get_results_path("MultiQC/logs/init.log")
    shell:
        r"""
        mkdir -p {output.multiqc_dir}
        mkdir -p {output.log_dir}
        echo "Created MultiQC directory structure" > {output.init_log}
        """

# Rule to generate MultiQC config file
rule generate_multiqc_config:
    input:
        init_log = get_results_path("MultiQC/logs/init.log"),
        custom_config = multiqc_config_file if use_custom_multiqc_config else []
    output:
        config = get_results_path("MultiQC/multiqc_config.yaml")
    log:
        get_results_path("MultiQC/logs/generate_multiqc_config.log")
    shell:
        r"""
        # Create logs directory
        mkdir -p $(dirname {log})
        
        echo "Generating MultiQC configuration file" > {log}
        
        if [ "{use_custom_multiqc_config}" = "True" ]; then
            echo "Using custom config file: {input.custom_config}" >> {log}
            cp {input.custom_config} {output.config}
        else
            echo "Creating default MultiQC config" >> {log}
            cat > {output.config} << 'EOF'
# MultiQC configuration
fn_ignore_files:
    - '*unique.txt.summary'
# Use first part of path as sample name
sample_names_regexes:
  - '^{config[results_dir]}/([^/]+)/'

# Clean up file extensions for better sample names

extra_fn_clean_exts:
  - '.*'
  - '_*'
  - '.1'
  - '.picard'
  - '_qualimap_bam'
  - '.1'
  - '_R1'
  - '_qualimap_rnaseq'  

fn_ignore_files:
    - '*_unique.txt.summary'

# Define specific patterns to search for
sp:
  fastp:
    fn: '*.fastp.*.json'
EOF
        fi
        
        echo "MultiQC configuration file generated successfully" >> {log}
        """

# Rule for paired-end samples
rule multiqc_paired_report:
    input:
        config_file = get_results_path("MultiQC/multiqc_config.yaml")
    output:
        html = get_results_path("MultiQC/multiqc_paired_report.html"),
        data_dir = directory(get_results_path("MultiQC/multiqc_paired_report_data"))
    params:
        qc_dirs = " ".join([os.path.join(config['results_dir'], sample, "qc") for sample in paired_samples]) if paired_samples else "",
        results_dir = get_results_path("MultiQC")
    log:
        get_results_path("MultiQC/logs/multiqc_paired.log")
    conda:
        "envs/multiqc.yaml"
    singularity:
        config.get("singularity_image", "docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0")
    shell:
        r"""
        # Create logs directory
        mkdir -p $(dirname {log})
        
        # Run MultiQC for paired-end samples
        echo "Running MultiQC on paired-end sample QC directories" > {log}
        
        if [ -n "{params.qc_dirs}" ]; then
            echo "Searching directories: {params.qc_dirs}" >> {log}
            multiqc {params.qc_dirs} \
                --ignore '*_unique.txt.summary' \
                --force \
                --config {input.config_file} \
                -o {params.results_dir} \
                -f \
                -n $(basename {output.html}) \
                --title "Paired-End Samples QC Report" \
                2>> {log}
            echo "MultiQC report for paired-end samples generated: {output.html}" >> {log}
        else
            echo "No paired-end samples found. Creating empty report." >> {log}
            mkdir -p {params.results_dir}
            echo "<html><body><h1>No paired-end samples found</h1></body></html>" > {output.html}
            mkdir -p {output.data_dir}
        fi
        """

# Rule for single-end samples
rule multiqc_single_report:
    input:
        config_file = get_results_path("MultiQC/multiqc_config.yaml")
    output:
        html = get_results_path("MultiQC/multiqc_single_report.html"),
        data_dir = directory(get_results_path("MultiQC/multiqc_single_report_data"))
    params:
        qc_dirs = " ".join([os.path.join(config['results_dir'], sample, "qc") for sample in single_samples]) if single_samples else "",
        results_dir = get_results_path("MultiQC")
    log:
        get_results_path("MultiQC/logs/multiqc_single.log")
    conda:
        "envs/multiqc.yaml"
    singularity:
        config.get("singularity_image", "docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0")
    shell:
        r"""
        # Create logs directory
        mkdir -p $(dirname {log})
        
        # Run MultiQC for single-end samples
        echo "Running MultiQC on single-end sample QC directories" > {log}
        
        if [ -n "{params.qc_dirs}" ]; then
            echo "Searching directories: {params.qc_dirs}" >> {log}
            multiqc {params.qc_dirs} \
                --ignore '*_unique.txt.summary' \
                --force \
                --config {input.config_file} \
                -o {params.results_dir} \
                -f \
                -n $(basename {output.html}) \
                --title "Single-End Samples QC Report" \
                2>> {log}
            echo "MultiQC report for single-end samples generated: {output.html}" >> {log}
        else
            echo "No single-end samples found. Creating empty report." >> {log}
            mkdir -p {params.results_dir}
            echo "<html><body><h1>No single-end samples found</h1></body></html>" > {output.html}
            mkdir -p {output.data_dir}
        fi
        """

# Combined rule to run both reports
rule multiqc_all_reports:
    input:
        paired_report = get_results_path("MultiQC/multiqc_paired_report.html"),
        single_report = get_results_path("MultiQC/multiqc_single_report.html")
    output:
        flag = get_results_path("MultiQC/multiqc_reports_complete")
    log:
        get_results_path("MultiQC/logs/multiqc_all_reports.log")
    shell:
        r"""
        # Create logs directory
        mkdir -p $(dirname {log})
        echo "All MultiQC reports generated successfully" | tee {log}
        echo "Paired-end report: {input.paired_report}" >> {log}
        echo "Single-end report: {input.single_report}" >> {log}
        touch {output.flag}
        """
