"""
Rules for alignment with bowtie2 for paired-end reads
Including read group information for Picard compatibility
"""
import os
import re

# Rule for aligning paired-end reads with bowtie2
rule align_paired_end:
    wildcard_constraints:
        sample = "|".join([re.escape(s) for s in PAIRED_LOCAL_SAMPLES + PAIRED_FTP_SAMPLES])
    input:
        r1 = get_cleaned_fastq_r1,
        r2 = get_cleaned_fastq_r2,
        fastp_flag = get_processing_path("{sample}/fastp_paired_complete.flag")
    output:
        bam = get_processing_path("{sample}/{sample}.bam"),
        bai = get_processing_path("{sample}/{sample}.bam.bai"),
        stats = get_processing_path("{sample}/qc/bowtie2/{sample}.bowtie2_paired_stats.txt"),
        flag = get_processing_path("{sample}/alignment_paired_complete.flag")
    params:
        genome_index = config["processing_genome_index"],  # Use the copied genome index
        # Read group parameters
        rg_id = lambda wildcards: wildcards.sample,
        rg_sm = lambda wildcards: wildcards.sample,
        rg_lb = lambda wildcards: f"{wildcards.sample}_lib",
        rg_pu = lambda wildcards: f"{wildcards.sample}_unit",
        rg_pl = "ILLUMINA"
    log:
        get_processing_path("{sample}/logs/alignment_paired.log")
    benchmark:
        get_processing_path("{sample}/benchmarks/align_paired.benchmark.txt")
    threads: 
        config["cores_align"]
    resources:
        mem_mb = 8000
    # Container and environment options
    conda:
        "../../envs/alignment.yaml"
    singularity:
        config.get("singularity_image", "")
    shell:
        """
        # Create necessary directories
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {output.stats})
        
        # Start logging
        echo "Aligning paired-end reads for {wildcards.sample}" > {log}
        echo "Using genome index: {params.genome_index}" >> {log}
        echo "Input files: {input.r1} and {input.r2}" >> {log}
        echo "Output BAM: {output.bam}" >> {log}
        echo "Including read group information:" >> {log}
        echo "  RG ID: {params.rg_id}" >> {log}
        echo "  RG SM: {params.rg_sm}" >> {log}
        echo "  RG LB: {params.rg_lb}" >> {log}
        echo "  RG PU: {params.rg_pu}" >> {log}
        echo "  RG PL: {params.rg_pl}" >> {log}
        
        # Check if genome index exists
        if [[ ! -f {params.genome_index}.1.bt2 && ! -f {params.genome_index}.1.bt2l ]]; then
            echo "ERROR: Genome index files not found at {params.genome_index}" >> {log}
            echo "Please check that the genome_index path in config.yaml is correct" >> {log}
            echo "and that the index files are accessible from the Singularity container." >> {log}
            exit 1
        fi
        
        # Check if input files exist and have content
        if [[ ! -s {input.r1} ]]; then
            echo "ERROR: Input file {input.r1} does not exist or is empty!" >> {log}
            exit 1
        fi
        
        if [[ ! -s {input.r2} ]]; then
            echo "ERROR: Input file {input.r2} does not exist or is empty!" >> {log}
            exit 1
        fi
        
        # Check for command availability
        echo "Checking for bowtie2 availability..." >> {log}
        if ! command -v bowtie2 &> /dev/null; then
            echo "ERROR: bowtie2 command not found!" >> {log}
            echo "Please make sure bowtie2 is installed and in PATH in your container." >> {log}
            exit 1
        fi
        
        echo "Checking for samtools availability..." >> {log}
        if ! command -v samtools &> /dev/null; then
            echo "ERROR: samtools command not found!" >> {log}
            echo "Please make sure samtools is installed and in PATH in your container." >> {log}
            exit 1
        fi
        
        # Verify bowtie2 and samtools versions
        echo "Bowtie2 version:" >> {log}
        bowtie2 --version | head -n 1 >> {log}
        echo "Samtools version:" >> {log}
        samtools --version | head -n 1 >> {log}
        
        # Check if output file already exists and has content
        if [[ -s {output.bam} && -s {output.bai} ]]; then
            echo "Alignment output already exists for {wildcards.sample}, skipping alignment" >> {log}
        else
            echo "Running bowtie2 on {wildcards.sample}" >> {log}
            
            # Run alignment with bowtie2 with read group information
            echo "Command: bowtie2 --very-sensitive-local -p {threads} -x {params.genome_index} \\
                     -1 {input.r1} -2 {input.r2} \\
                     --rg-id '{params.rg_id}' \\
                     --rg 'SM:{params.rg_sm}' \\
                     --rg 'LB:{params.rg_lb}' \\
                     --rg 'PU:{params.rg_pu}' \\
                     --rg 'PL:{params.rg_pl}'" >> {log}
            
            # Run with output captured to QC file
            (bowtie2 --very-sensitive-local -p {threads} -x {params.genome_index} \\
                     -1 {input.r1} -2 {input.r2} \\
                     --rg-id '{params.rg_id}' \\
                     --rg 'SM:{params.rg_sm}' \\
                     --rg 'LB:{params.rg_lb}' \\
                     --rg 'PU:{params.rg_pu}' \\
                     --rg 'PL:{params.rg_pl}' \\
                     2> {output.stats} | \\
             samtools view -bSu -@ {threads} | \\
             samtools sort -@ {threads} -o {output.bam}) >> {log} 2>&1
             
            # Check if BAM file was created successfully
            if [[ ! -s {output.bam} ]]; then
                echo "ERROR: Failed to create BAM file. Check logs for details." >> {log}
                cat {output.stats} >> {log}
                exit 1
            fi
            
            # Create index for the BAM file
            echo "Indexing BAM file with {threads} threads..." >> {log}
            samtools index -@ {threads} {output.bam}
            
            echo "Alignment complete for {wildcards.sample}" >> {log}
        fi
        
        # Create flag file to indicate completion
        echo "Alignment complete for paired-end sample {wildcards.sample}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "Files created:" >> {output.flag}
        echo "- BAM file: {output.bam}" >> {output.flag}
        echo "- BAM index: {output.bai}" >> {output.flag}
        echo "- Bowtie2 stats: {output.stats}" >> {output.flag}
        """