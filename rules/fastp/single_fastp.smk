"""
Rules for quality filtering with fastp for single-end reads
"""
import os

# Rule for fastp processing of single-end reads
rule fastp_single_end:
    input:
        r = get_fastq_single,
        checksums = get_processing_path("{sample}/checksums_single_verified.flag")
    output:
        r = get_processing_path("{sample}/{sample}.cleaned.fastq.gz"),
        html = get_processing_path("{sample}/qc/fastp/{sample}.fastp.single.html"),
        json = get_processing_path("{sample}/qc/fastp/{sample}.fastp.single.json"),
        flag = get_processing_path("{sample}/fastp_single_complete.flag")
    log:
        get_processing_path("{sample}/logs/fastp_single.log")
    benchmark:
        get_processing_path("{sample}/benchmarks/fastp_single.benchmark.txt")
    threads: 
        config["cores_fastp"]
    resources:
        mem_mb = lambda wildcards: config.get("mem_fastp", 4000)
    # Container and environment options
    conda:
        "../../envs/fastp.yaml"
    singularity:
        config.get("singularity_image", "")
    shell:
        """
        # Create output directories
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.html})
        mkdir -p $(dirname {output.r})
        
        echo "Processing single-end fastq file for {wildcards.sample}" > {log}
        
        # Check if output files already exist and have content
        if [[ -s {output.r} && -s {output.html} && -s {output.json} ]]; then
            echo "Fastp outputs already exist for {wildcards.sample}, skipping processing" >> {log}
        else
            echo "Running fastp on {wildcards.sample}" >> {log}
            echo "Input file: {input.r}" >> {log}
            echo "Output file: {output.r}" >> {log}
            
            # Process files
            fastp -i {input.r} \\
                -o {output.r} \\
                -h {output.html} -j {output.json} \\
                -w {threads} \\
                2>> {log}
                
            # Check if output file was created successfully
            if [[ ! -s {output.r} ]]; then
                echo "ERROR: Failed to create output file" >> {log}
                exit 1
            fi
        fi
        
        # Create flag file to indicate completion
        echo "Fastp processing complete for single-end sample {wildcards.sample}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "Files processed:" >> {output.flag}
        echo "- Input: {input.r}" >> {output.flag}
        echo "- Output: {output.r}" >> {output.flag}
        echo "- HTML report: {output.html}" >> {output.flag}
        echo "- JSON report: {output.json}" >> {output.flag}
        """