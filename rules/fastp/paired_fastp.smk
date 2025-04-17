"""
Rules for quality filtering with fastp for paired-end reads
"""
import os

# Rule for fastp processing of paired-end reads
rule fastp_paired_end:
    input:
        r1 = get_fastq_r1,
        r2 = get_fastq_r2,
        checksums = get_processing_path("{sample}/checksums_paired_verified.flag")
    output:
        r1 = get_processing_path("{sample}/{sample}.1.cleaned.fastq.gz"),
        r2 = get_processing_path("{sample}/{sample}.2.cleaned.fastq.gz"),
        html = get_processing_path("{sample}/qc/fastp/{sample}.fastp.paired.html"),
        json = get_processing_path("{sample}/qc/fastp/{sample}.fastp.paired.json"),
        flag = get_processing_path("{sample}/fastp_paired_complete.flag")
    log:
        get_processing_path("{sample}/logs/fastp_paired.log")
    benchmark:
        get_processing_path("{sample}/benchmarks/fastp_paired.benchmark.txt")
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
        mkdir -p $(dirname {output.r1})
        
        echo "Processing paired-end fastq files for {wildcards.sample}" > {log}
        
        # Check if output files already exist and have content
        if [[ -s {output.r1} && -s {output.r2} && -s {output.html} && -s {output.json} ]]; then
            echo "Fastp outputs already exist for {wildcards.sample}, skipping processing" >> {log}
        else
            echo "Running fastp on {wildcards.sample}" >> {log}
            echo "Input files: {input.r1} and {input.r2}" >> {log}
            echo "Output files: {output.r1} and {output.r2}" >> {log}
            
            # Process files
            fastp -i {input.r1} -I {input.r2} \\
                -o {output.r1} -O {output.r2} \\
                -h {output.html} -j {output.json} \\
                -w {threads} \\
                2>> {log}
                
            # Check if output files were created successfully
            if [[ ! -s {output.r1} || ! -s {output.r2} ]]; then
                echo "ERROR: Failed to create output files" >> {log}
                exit 1
            fi
        fi
        
        # Create flag file to indicate completion
        echo "Fastp processing complete for paired-end sample {wildcards.sample}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "Files processed:" >> {output.flag}
        echo "- Input R1: {input.r1}" >> {output.flag}
        echo "- Input R2: {input.r2}" >> {output.flag}
        echo "- Output R1: {output.r1}" >> {output.flag}
        echo "- Output R2: {output.r2}" >> {output.flag}
        echo "- HTML report: {output.html}" >> {output.flag}
        echo "- JSON report: {output.json}" >> {output.flag}
        """