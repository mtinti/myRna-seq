"""
Rules for basic QC analysis of aligned and processed BAM files
Includes:
- samtools flagstat
- samtools stats
- qualimap bamqc
"""
import os

# Rule for running samtools flagstat
rule samtools_flagstat:
    input:
        bam = get_picard_bam,
        bai = lambda wildcards: f"{get_picard_bam(wildcards)}.bai",
        markduplicates_flag = get_processing_path("{sample}/markduplicates_complete.flag")
    output:
        flagstat = get_processing_path("{sample}/qc/flagstat/{sample}.flagstat.txt")
    log:
        get_processing_path("{sample}/logs/flagstat.log")
    benchmark:
        get_processing_path("{sample}/benchmarks/flagstat.benchmark.txt")
    threads: 
        config["cores_flagstat"]
    conda:
        "../../envs/samtools.yaml"
    singularity:
        config.get("singularity_image", "")
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.flagstat})
        mkdir -p $(dirname {log})
        
        # Log start
        echo "Running samtools flagstat on {wildcards.sample}" > {log}
        echo "Input BAM: {input.bam}" >> {log}
        echo "Output: {output.flagstat}" >> {log}
        
        # Run flagstat
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat} 2>> {log}
        
        # Check if output was created successfully
        if [[ ! -s {output.flagstat} ]]; then
            echo "ERROR: Failed to create flagstat output" >> {log}
            exit 1
        fi
        
        echo "Flagstat completed successfully" >> {log}
        """

# Rule for running samtools stats
rule samtools_stats:
    input:
        bam = get_picard_bam,
        bai = lambda wildcards: f"{get_picard_bam(wildcards)}.bai",
        markduplicates_flag = get_processing_path("{sample}/markduplicates_complete.flag")
    output:
        stats = get_processing_path("{sample}/qc/stats/{sample}.stats.txt")
    log:
        get_processing_path("{sample}/logs/stats.log")
    benchmark:
        get_processing_path("{sample}/benchmarks/stats.benchmark.txt")
    threads: 
        config["cores_flagstat"]
    conda:
        "../../envs/samtools.yaml"
    singularity:
        config.get("singularity_image", "")
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.stats})
        mkdir -p $(dirname {log})
        
        # Log start
        echo "Running samtools stats on {wildcards.sample}" > {log}
        echo "Input BAM: {input.bam}" >> {log}
        echo "Output: {output.stats}" >> {log}
        
        # Run stats
        samtools stats -@ {threads} {input.bam} > {output.stats} 2>> {log}
        
        # Check if output was created successfully
        if [[ ! -s {output.stats} ]]; then
            echo "ERROR: Failed to create stats output" >> {log}
            exit 1
        fi
        
        echo "Stats completed successfully" >> {log}
        """

# Rule for running qualimap bamqc
rule qualimap_bamqc:
    input:
        bam = get_picard_bam,
        bai = lambda wildcards: f"{get_picard_bam(wildcards)}.bai",
        markduplicates_flag = get_processing_path("{sample}/markduplicates_complete.flag")
    output:
        dir = directory(get_processing_path("{sample}/qc/{sample}_qualimap_bam/")),
        flag = get_processing_path("{sample}/bamqc_complete.flag")
    log:
        stderr = get_processing_path("{sample}/logs/qualimap_bamqc.log"),
        stdout = get_processing_path("{sample}/logs/qualimap_bamqc.stdout.log")
    benchmark:
        get_processing_path("{sample}/benchmarks/qualimap_bamqc.benchmark.txt")
    params:
        memory = config.get("qualimap_memory", "10G")
    threads: 
        config["cores_flagstat"]
    conda:
        "../../envs/qualimap.yaml"
    singularity:
        config.get("singularity_image", "")
    shell:
        """
        # Create output directory
        mkdir -p {output.dir}
        mkdir -p $(dirname {log.stderr})
        
        # Log start
        echo "Running qualimap bamqc on {wildcards.sample}" > {log.stderr}
        echo "Input BAM: {input.bam}" >> {log.stderr}
        echo "Output directory: {output.dir}" >> {log.stderr}
        echo "Java memory: {params.memory}" >> {log.stderr}
        
        # Run qualimap bamqc
        qualimap bamqc --java-mem-size={params.memory} \\
            -bam {input.bam} \\
            -outdir {output.dir} \\
            -outformat HTML \\
            -nt {threads} > {log.stdout} 2>> {log.stderr}
        
        # Check for successful completion by looking for report.html
        if [[ ! -f {output.dir}/qualimapReport.html ]]; then
            echo "ERROR: qualimap bamqc did not produce report.html" >> {log.stderr}
            exit 1
        fi
        
        # Create flag file to indicate completion
        echo "BAM QC complete for sample {wildcards.sample}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "Files generated:" >> {output.flag}
        echo "- Qualimap BAM QC report: {output.dir}/qualimapReport.html" >> {output.flag}
        
        echo "Qualimap bamqc completed successfully" >> {log.stderr}
        """