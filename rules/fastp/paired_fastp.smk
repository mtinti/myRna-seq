"""
Rules for quality filtering with fastp for paired-end and split interleaved reads
"""
import os

def get_fastp_input(wildcards):
    """Return input files for the fastp rule.

    Both ``split_complete`` and ``checksums`` keys are always present so the
    shell block can reference them without conditional logic failing when mixing
    sample types.
    """

    if SAMPLES_DF.loc[wildcards.sample, 'read_type'] == 'interleaved':
        return {
            'r1': get_processing_path(f"{wildcards.sample}/{wildcards.sample}.1.fastq.gz"),
            'r2': get_processing_path(f"{wildcards.sample}/{wildcards.sample}.2.fastq.gz"),
            'split_complete': get_processing_path(f"{wildcards.sample}/split_interleaved_complete.flag"),
            'checksums': ''
        }
    else:
        return {
            'r1': get_processing_path(f"{wildcards.sample}/{wildcards.sample}.1.fastq.gz"),
            'r2': get_processing_path(f"{wildcards.sample}/{wildcards.sample}.2.fastq.gz"),
            'split_complete': '',
            'checksums': get_processing_path(f"{wildcards.sample}/checksums_paired_verified.flag")
        }

# Rule for fastp processing of paired-end and split interleaved reads
rule fastp_paired_end:
    input:
        unpack(get_fastp_input)
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
        
        # Determine input file type
        if [[ -n "{input.split_complete}" && -f "{input.split_complete}" ]]; then
            echo "Processing split interleaved files for {wildcards.sample}" > {log}
            input_r1={input.r1}
            input_r2={input.r2}
        else
            echo "Processing paired-end fastq files for {wildcards.sample}" > {log}
            input_r1={input.r1}
            input_r2={input.r2}
        fi
        
        # Check if output files already exist and have content
        if [[ -s {output.r1} && -s {output.r2} && -s {output.html} && -s {output.json} ]]; then
            echo "Fastp outputs already exist for {wildcards.sample}, skipping processing" >> {log}
        else
            echo "Running fastp on {wildcards.sample}" >> {log}
            echo "Input files: $input_r1 and $input_r2" >> {log}
            echo "Output files: {output.r1} and {output.r2}" >> {log}
            
            # Process files
            fastp -i $input_r1 -I $input_r2 \\
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
        if [[ -n "{input.split_complete}" && -f "{input.split_complete}" ]]; then
            sample_type="split interleaved"
        else
            sample_type="paired-end"
        fi
        echo "Fastp processing complete for $sample_type sample {wildcards.sample}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "Files processed:" >> {output.flag}
        echo "- Input R1: $input_r1" >> {output.flag}
        echo "- Input R2: $input_r2" >> {output.flag}
        echo "- Output R1: {output.r1}" >> {output.flag}
        echo "- Output R2: {output.r2}" >> {output.flag}
        echo "- HTML report: {output.html}" >> {output.flag}
        echo "- JSON report: {output.json}" >> {output.flag}
        """