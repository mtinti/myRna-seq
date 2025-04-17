"""
Rules for counting reads in genomic features for paired-end samples using featureCounts
Handles both all reads and uniquely mapped reads
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

# Rule for counting all reads in paired-end data
rule featurecounts_paired_all:
    wildcard_constraints:
        run_tag = "|".join([re.escape(rt) for rt in PAIRED_RUN_TAGS]) if PAIRED_RUN_TAGS else "^$"
    input:
        bam = get_picard_bam,
        bai = lambda wildcards: f"{get_picard_bam(wildcards)}.bai",
        coverage_complete = get_processing_path("{run_tag}/coverage_paired_complete.flag")
    output:
        counts = get_processing_path("{run_tag}/{run_tag}_counts_paired_all.txt"),
        summary = get_processing_path("{run_tag}/qc/feature_counts/{run_tag}_counts_paired_all.txt.summary")
    params:
        gtf = config["processing_gtf_file"],
        feature_type = config["feature_type"],
        extra = "-p -B -C -M -O --countReadPairs"  # Paired-end specific options
    log:
        get_processing_path("{run_tag}/logs/featurecounts_paired_all.log")
    benchmark:
        get_processing_path("{run_tag}/benchmarks/featurecounts_paired_all.benchmark.txt")
    threads: 
        config["cores_featurecounts"]
    conda:
        "../../envs/subread.yaml"
    singularity:
        config.get("singularity_image", "")
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.counts})
        mkdir -p $(dirname {output.summary})
        mkdir -p $(dirname {log})
        
        # Log start
        echo "Running featureCounts (paired-end, all reads) for {wildcards.run_tag}" > {log}
        echo "Input BAM: {input.bam}" >> {log}
        echo "Output: {output.counts}" >> {log}
        echo "Summary: {output.summary}" >> {log}
        echo "GTF file: {params.gtf}" >> {log}
        echo "Feature type: {params.feature_type}" >> {log}
        echo "Extra options: {params.extra}" >> {log}
        
        # Check if GTF file exists
        if [[ ! -f {params.gtf} ]]; then
            echo "ERROR: GTF file {params.gtf} not found" >> {log}
            exit 1
        fi
        
        # Run featureCounts for paired-end data, all reads
        featureCounts {params.extra} \\
            -T {threads} \\
            -t {params.feature_type} \\
            -g gene_id \\
            -a {params.gtf} \\
            -o {output.counts} \\
            {input.bam} \\
            2>> {log}
        
        # Check if output was created successfully
        if [[ ! -s {output.counts} ]]; then
            echo "ERROR: Failed to create counts file" >> {log}
            exit 1
        fi
        
        # featureCounts automatically creates a summary file with .summary extension
        auto_summary="{output.counts}.summary"
        if [[ -f $auto_summary ]]; then
            echo "Moving summary file to QC directory" >> {log}
            mv $auto_summary {output.summary}
        else
            echo "ERROR: Summary file not found at $auto_summary" >> {log}
            exit 1
        fi
        
        echo "featureCounts (paired-end, all reads) completed successfully" >> {log}
        """

# Rule for counting uniquely mapped reads in paired-end data
rule featurecounts_paired_unique:
    wildcard_constraints:
        run_tag = "|".join([re.escape(rt) for rt in PAIRED_RUN_TAGS]) if PAIRED_RUN_TAGS else "^$"
    input:
        bam = get_picard_bam,
        bai = lambda wildcards: f"{get_picard_bam(wildcards)}.bai",
        coverage_complete = get_processing_path("{run_tag}/coverage_paired_complete.flag"),
        counts_all = get_processing_path("{run_tag}/{run_tag}_counts_paired_all.txt")
    output:
        counts = get_processing_path("{run_tag}/{run_tag}_counts_paired_unique.txt"),
        summary = get_processing_path("{run_tag}/qc/feature_counts/{run_tag}_counts_paired_unique.txt.summary"),
        flag = get_processing_path("{run_tag}/featurecounts_paired_complete.flag")
    params:
        gtf = config["processing_gtf_file"],
        feature_type = config["feature_type"], 
        extra = "-p -B -C -O --countReadPairs",  # Paired-end specific options
        min_mapping_quality = config.get("min_mapping_quality", 2)
    log:
        get_processing_path("{run_tag}/logs/featurecounts_paired_unique.log")
    benchmark:
        get_processing_path("{run_tag}/benchmarks/featurecounts_paired_unique.benchmark.txt")
    threads: 
        config["cores_featurecounts"]
    conda:
        "../../envs/subread.yaml"
    singularity:
        config.get("singularity_image", "")
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.counts})
        mkdir -p $(dirname {output.summary})
        mkdir -p $(dirname {log})
        
        # Log start
        echo "Running featureCounts (paired-end, unique reads) for {wildcards.run_tag}" > {log}
        echo "Input BAM: {input.bam}" >> {log}
        echo "Output: {output.counts}" >> {log}
        echo "Summary: {output.summary}" >> {log}
        echo "GTF file: {params.gtf}" >> {log}
        echo "Feature type: {params.feature_type}" >> {log}
        echo "Extra options: {params.extra}" >> {log}
        echo "Min mapping quality: {params.min_mapping_quality}" >> {log}
        
        # Check if GTF file exists
        if [[ ! -f {params.gtf} ]]; then
            echo "ERROR: GTF file {params.gtf} not found" >> {log}
            exit 1
        fi
        
        # Run featureCounts for paired-end data, uniquely mapped reads
        featureCounts {params.extra} \\
            -T {threads} \\
            -t {params.feature_type} \\
            -g gene_id \\
            -Q {params.min_mapping_quality} \\
            -a {params.gtf} \\
            -o {output.counts} \\
            {input.bam} \\
            2>> {log}
        
        # Check if output was created successfully
        if [[ ! -s {output.counts} ]]; then
            echo "ERROR: Failed to create counts file" >> {log}
            exit 1
        fi
        
        # featureCounts automatically creates a summary file with .summary extension
        auto_summary="{output.counts}.summary"
        if [[ -f $auto_summary ]]; then
            echo "Moving summary file to QC directory" >> {log}
            mv $auto_summary {output.summary}
        else
            echo "ERROR: Summary file not found at $auto_summary" >> {log}
            exit 1
        fi
        
        # Create flag file to indicate completion
        echo "Feature counting completed for paired-end sample {wildcards.run_tag}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "---------------------------------" >> {output.flag}
        echo "All reads counts: {input.counts_all}" >> {output.flag}
        echo "Unique reads counts: {output.counts}" >> {output.flag}
        
        echo "featureCounts (paired-end, unique reads) completed successfully" >> {log}
        """