"""
Rule for marking duplicates with Picard - modified to use merged BAM files
Creates .picard.bam files with duplicates either marked or removed based on config
"""
import os

# Unified rule for marking duplicates using merged BAM files
rule mark_duplicates:
    input:
        # Use merged BAM files as input
        bam = get_processing_path("merged/{run_tag}/{run_tag}.bam"),
        bai = get_processing_path("merged/{run_tag}/{run_tag}.bam.bai"),
        merge_flag = get_processing_path("merged/{run_tag}/merge_complete.flag")
    output:
        bam = get_processing_path("{run_tag}/{run_tag}.picard.bam"),
        bai = get_processing_path("{run_tag}/{run_tag}.picard.bam.bai"),
        metrics = get_processing_path("{run_tag}/qc/markduplicates/{run_tag}.markduplicates_metrics.txt"),
        flag = get_processing_path("{run_tag}/markduplicates_complete.flag")
    params:
        remove_duplicates = lambda wildcards: "true" if str(config.get("remove_duplicates", False)).lower() == "true" else "false"
    log:
        get_processing_path("{run_tag}/logs/markduplicates.log")
    benchmark:
        get_processing_path("{run_tag}/benchmarks/markduplicates.benchmark.txt")
    threads: 
        config["cores_default"]
    resources:
        mem_mb = 8000
    # Container and environment options
    conda:
        "../../envs/picard.yaml"
    singularity:
        config.get("singularity_image", "")
    shell:
        """
        # Create output directories
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {output.metrics})
        mkdir -p $(dirname {log})
        
        # Log start
        echo "Running MarkDuplicates for {wildcards.run_tag}" > {log}
        echo "Input BAM: {input.bam}" >> {log}
        echo "Output BAM: {output.bam}" >> {log}
        echo "Metrics file: {output.metrics}" >> {log}
        echo "Remove duplicates: {params.remove_duplicates}" >> {log}
        echo "Using {threads} threads for processing" >> {log}
        
        # Check if output files already exist
        if [[ -s {output.bam} && -s {output.bai} && -s {output.metrics} ]]; then
            echo "MarkDuplicates outputs already exist for {wildcards.run_tag}, skipping processing" >> {log}
            exit 0
        fi
        
        # Check for picard availability
        echo "Checking for picard availability..." >> {log}
        if ! command -v picard &> /dev/null; then
            if command -v picard.jar &> /dev/null; then
                PICARD_CMD="picard.jar"
            elif [ -f $PICARD_JAR ]; then
                PICARD_CMD="java -jar $PICARD_JAR"
            else
                echo "ERROR: picard command not found" >> {log}
                exit 1
            fi
        else
            PICARD_CMD="picard"
        fi
        
        echo "Using picard command: $PICARD_CMD" >> {log}
        
        # Run MarkDuplicates
        $PICARD_CMD MarkDuplicates \\
            I={input.bam} \\
            O={output.bam} \\
            M={output.metrics} \\
            REMOVE_DUPLICATES={params.remove_duplicates} \\
            ASSUME_SORT_ORDER=coordinate \\
            VALIDATION_STRINGENCY=LENIENT \\
            2>> {log}
        
        # Check if BAM file was created successfully
        if [[ ! -s {output.bam} ]]; then
            echo "ERROR: Failed to create BAM file" >> {log}
            exit 1
        fi
        
        # Index the BAM file with threads
        echo "Indexing BAM file with {threads} threads..." >> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        
        # Check if index was created successfully
        if [[ ! -s {output.bai} ]]; then
            echo "ERROR: Failed to create BAM index" >> {log}
            exit 1
        fi
        
        # Create flag file to indicate completion
        echo "MarkDuplicates complete for run tag {wildcards.run_tag}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "Files created:" >> {output.flag}
        echo "- BAM file: {output.bam}" >> {output.flag}
        echo "- BAM index: {output.bai}" >> {output.flag}
        echo "- Metrics file: {output.metrics}" >> {output.flag}
        echo "Remove duplicates: {params.remove_duplicates}" >> {output.flag}
        
        echo "MarkDuplicates complete for {wildcards.run_tag}" >> {log}
        """