"""
Rules for verifying file checksums for paired-end samples
Fails if checksums don't match (can be overridden with --keep-going)
Uses shared checksum verification functions
"""
import os

# Rule for checking MD5 checksums for paired-end data
rule verify_checksums_paired:
    input:
        r1 = get_fastq_r1,
        r2 = get_fastq_r2,
        acquisition_flag = get_processing_path("{sample}/acquisition_paired_complete.flag"),
        functions = "rules/checksum/checksum_functions.sh"
    output:
        flag = get_processing_path("{sample}/checksums_paired_verified.flag")
    params:
        checksum_r1 = lambda wildcards: SAMPLES_DF.loc[wildcards.sample, 'checksum_1'] if 'checksum_1' in SAMPLES_DF.columns else "",
        checksum_r2 = lambda wildcards: SAMPLES_DF.loc[wildcards.sample, 'checksum_2'] if 'checksum_2' in SAMPLES_DF.columns else ""
    log:
        get_processing_path("{sample}/logs/checksum_paired.log")
    benchmark:
        get_processing_path("{sample}/benchmarks/checksum_paired.benchmark.txt")
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
        echo "Verifying checksums for paired-end files for sample {wildcards.sample}" > {log}
        echo "R1: {input.r1}" >> {log}
        echo "R2: {input.r2}" >> {log}
        
        # Initialize flag file with timestamp and sample info
        echo "CHECKSUM VERIFICATION RESULTS FOR PAIRED-END SAMPLE {wildcards.sample}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "----------------------------------------" >> {output.flag}
        
        # Source the common functions
        source {input.functions}
        
        # Check R1 file
        check_md5sum "{input.r1}" "{params.checksum_r1}" "R1" "{output.flag}" "{log}"
        
        # Check R2 file
        check_md5sum "{input.r2}" "{params.checksum_r2}" "R2" "{output.flag}" "{log}"
        
        # Mark completion
        echo "----------------------------------------" >> {output.flag}
        echo "STATUS: COMPLETE" >> {output.flag}
        echo "All checksums verified successfully" >> {output.flag}
        
        echo "Checksum verification completed successfully" >> {log}
        """