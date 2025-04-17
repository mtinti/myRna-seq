"""
Rules for verifying file checksums for single-end samples
Fails if checksums don't match (can be overridden with --keep-going)
Uses shared checksum verification functions
"""
import os

# Rule for checking MD5 checksums for single-end data
rule verify_checksums_single:
    input:
        r = get_fastq_single,
        acquisition_flag = get_processing_path("{sample}/acquisition_single_complete.flag"),
        functions = "rules/checksum/checksum_functions.sh"
    output:
        flag = get_processing_path("{sample}/checksums_single_verified.flag")
    params:
        checksum = lambda wildcards: SAMPLES_DF.loc[wildcards.sample, 'checksum_1'] if 'checksum_1' in SAMPLES_DF.columns else ""
    log:
        get_processing_path("{sample}/logs/checksum_single.log")
    benchmark:
        get_processing_path("{sample}/benchmarks/checksum_single.benchmark.txt")
    resources:
        mem_mb = lambda wildcards: config.get("mem_checksum", 1000)
    conda:
        "../../envs/core.yaml"
    singularity:
        config.get("singularity_image", "")
    shell:
        """
        # Create output directories
        mkdir -p $(dirname {log})
        
        # Start logging
        echo "Verifying checksum for single-end file for sample {wildcards.sample}" > {log}
        echo "File: {input.r}" >> {log}
        
        # Initialize flag file with timestamp and sample info
        echo "CHECKSUM VERIFICATION RESULTS FOR SINGLE-END SAMPLE {wildcards.sample}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "----------------------------------------" >> {output.flag}
        
        # Source the common functions
        source {input.functions}
        
        # Check the single file
        check_md5sum "{input.r}" "{params.checksum}" "R1" "{output.flag}" "{log}"
        
        # Mark completion
        echo "----------------------------------------" >> {output.flag}
        echo "STATUS: COMPLETE" >> {output.flag}
        echo "All checksums verified successfully" >> {output.flag}
        
        echo "Checksum verification completed successfully" >> {log}
        """