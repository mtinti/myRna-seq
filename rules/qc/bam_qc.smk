"""
Rules for basic QC analysis of aligned and processed BAM files
Includes:
- samtools flagstat
- samtools stats
- qualimap bamqc
"""
import os
import re

# List of run tags that are not nanopore
NON_NANOPORE_RUN_TAGS = [rt for rt in EFFECTIVE_SAMPLES_TO_PROCESS if not is_nanopore(rt)]

# Rule for running samtools flagstat
rule samtools_flagstat:
    input:
        bam = get_picard_bam,
        bai = lambda wildcards: f"{get_picard_bam(wildcards)}.bai",
        markduplicates_flag = get_markduplicates_flag
    output:
        flagstat = get_processing_path("{run_tag}/qc/flagstat/{run_tag}.flagstat.txt")
    log:
        get_processing_path("{run_tag}/logs/flagstat.log")
    benchmark:
        get_processing_path("{run_tag}/benchmarks/flagstat.benchmark.txt")
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
        echo "Running samtools flagstat on {wildcards.run_tag}" > {log}
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
        markduplicates_flag = get_markduplicates_flag
    output:
        stats = get_processing_path("{run_tag}/qc/stats/{run_tag}.stats.txt")
    log:
        get_processing_path("{run_tag}/logs/stats.log")
    benchmark:
        get_processing_path("{run_tag}/benchmarks/stats.benchmark.txt")
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
        echo "Running samtools stats on {wildcards.run_tag}" > {log}
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
    wildcard_constraints:
        run_tag = "|".join([re.escape(rt) for rt in NON_NANOPORE_RUN_TAGS]) if NON_NANOPORE_RUN_TAGS else "^$"
    input:
        bam = get_picard_bam,
        bai = lambda wildcards: f"{get_picard_bam(wildcards)}.bai",
        markduplicates_flag = get_markduplicates_flag
    output:
        dir = directory(get_processing_path("{run_tag}/qc/{run_tag}_qualimap_bam/")),
        flag = get_processing_path("{run_tag}/bamqc_complete.flag")
    log:
        stderr = get_processing_path("{run_tag}/logs/qualimap_bamqc.log"),
        stdout = get_processing_path("{run_tag}/logs/qualimap_bamqc.stdout.log")
    benchmark:
        get_processing_path("{run_tag}/benchmarks/qualimap_bamqc.benchmark.txt")
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
        echo "Running qualimap bamqc on {wildcards.run_tag}" > {log.stderr}
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
        echo "BAM QC complete for run tag {wildcards.run_tag}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "Files generated:" >> {output.flag}
        echo "- Qualimap BAM QC report: {output.dir}/qualimapReport.html" >> {output.flag}
        
        echo "Qualimap bamqc completed successfully" >> {log.stderr}
        """