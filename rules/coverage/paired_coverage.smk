"""
Rules for generating coverage tracks for paired-end samples using bamCoverage from deepTools
Creates both all-read and unique-read coverage tracks
Updated to work with run tags
"""
import os
import re

# Helper function to check if a run tag is for paired-end samples
def is_valid_paired_run_tag(run_tag):
    # For standard samples that are their own run tag
    if run_tag in SAMPLES_DF.index and run_tag not in RUN_TAGS:
        return is_paired_end(run_tag)
        
    # For actual run tags that group multiple samples
    if run_tag in RUN_TAG_SAMPLES:
        samples = RUN_TAG_SAMPLES[run_tag]
        if samples:
            # Check the first sample in the group to determine type
            return is_paired_end(samples[0])
    
    # Default case - assume it's not valid for this rule
    return False

# Define which run tags are paired-end
PAIRED_RUN_TAGS = [rt for rt in EFFECTIVE_SAMPLES_TO_PROCESS if is_valid_paired_run_tag(rt)]

rule generate_all_reads_coverage_paired:
    wildcard_constraints:
        run_tag = "|".join([re.escape(rt) for rt in PAIRED_RUN_TAGS]) if PAIRED_RUN_TAGS else "^$"
    input:
        bam = get_picard_bam,
        bai = lambda wildcards: f"{get_picard_bam(wildcards)}.bai",
        qc_complete = get_processing_path("{run_tag}/qc_paired_complete.flag")
    output:
        bigwig = get_processing_path("{run_tag}/{run_tag}_all.bw")
    log:
        get_processing_path("{run_tag}/logs/coverage_all.log")
    benchmark:
        get_processing_path("{run_tag}/benchmarks/coverage_all.benchmark.txt")
    params:
        bin_size = config.get("coverage_bin_size", 50),
        normalize = config.get("coverage_normalize", "RPKM")
    threads: 
        config["cores_coverage"]
    conda:
        "../../envs/deeptools.yaml"
    singularity:
        config.get("singularity_image", "")
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.bigwig})
        mkdir -p $(dirname {log})
        
        # Log start
        echo "Generating coverage track (all reads) for paired-end sample {wildcards.run_tag}" > {log}
        echo "Input BAM: {input.bam}" >> {log}
        echo "Output: {output.bigwig}" >> {log}
        echo "Bin size: {params.bin_size}" >> {log}
        echo "Normalization: {params.normalize}" >> {log}
        
        # Generate coverage track for all reads
        bamCoverage --bam {input.bam} \\
            --outFileName {output.bigwig} \\
            --binSize {params.bin_size} \\
            --normalizeUsing {params.normalize} \\
            --numberOfProcessors {threads} \\
            --ignoreDuplicates \\
            >> {log} 2>&1
        
        # Check if output was created successfully
        if [[ ! -s {output.bigwig} ]]; then
            echo "ERROR: Failed to create coverage track" >> {log}
            exit 1
        fi
        
        echo "Coverage track generation complete" >> {log}
        """

# Rule for generating coverage track with uniquely mapped reads only for paired-end samples
rule generate_unique_reads_coverage_paired:
    wildcard_constraints:
        run_tag = "|".join([re.escape(rt) for rt in PAIRED_RUN_TAGS]) if PAIRED_RUN_TAGS else "^$"
    input:
        bam = get_picard_bam,
        bai = lambda wildcards: f"{get_picard_bam(wildcards)}.bai",
        qc_complete = get_processing_path("{run_tag}/qc_paired_complete.flag"),
        all_bigwig = get_processing_path("{run_tag}/{run_tag}_all.bw")
    output:
        bigwig = get_processing_path("{run_tag}/{run_tag}_unique.bw"),
        flag = get_processing_path("{run_tag}/coverage_paired_complete.flag")
    log:
        get_processing_path("{run_tag}/logs/coverage_unique.log")
    benchmark:
        get_processing_path("{run_tag}/benchmarks/coverage_unique.benchmark.txt")
    params:
        bin_size = config.get("coverage_bin_size", 50),
        normalize = config.get("coverage_normalize", "RPKM"),
        min_mapping_quality = config.get("min_mapping_quality", 2)
    threads: 
        config["cores_coverage"]
    conda:
        "../../envs/deeptools.yaml"
    singularity:
        config.get("singularity_image", "")
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.bigwig})
        mkdir -p $(dirname {log})
        
        # Log start
        echo "Generating coverage track (unique reads) for paired-end sample {wildcards.run_tag}" > {log}
        echo "Input BAM: {input.bam}" >> {log}
        echo "Output: {output.bigwig}" >> {log}
        echo "Bin size: {params.bin_size}" >> {log}
        echo "Normalization: {params.normalize}" >> {log}
        echo "Min mapping quality: {params.min_mapping_quality}" >> {log}
        
        # Generate coverage track for uniquely mapped reads
        bamCoverage --bam {input.bam} \\
            --outFileName {output.bigwig} \\
            --binSize {params.bin_size} \\
            --normalizeUsing {params.normalize} \\
            --minMappingQuality {params.min_mapping_quality} \\
            --numberOfProcessors {threads} \\
            --ignoreDuplicates \\
            >> {log} 2>&1
        
        # Check if output was created successfully
        if [[ ! -s {output.bigwig} ]]; then
            echo "ERROR: Failed to create coverage track" >> {log}
            exit 1
        fi
        
        # Create flag file to indicate completion
        echo "Coverage track generation completed for paired-end sample {wildcards.run_tag}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "---------------------------------" >> {output.flag}
        echo "All reads track: {input.all_bigwig}" >> {output.flag}
        echo "Unique reads track: {output.bigwig}" >> {output.flag}
        
        echo "Coverage track generation complete" >> {log}
        """