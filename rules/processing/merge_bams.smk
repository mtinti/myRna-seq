"""
Rules for merging BAM files from samples with the same run tag
Handles both actual merging and simple copying for standalone samples
"""
import os

# Rule for merging BAM files from samples with the same run tag (or copying for single-sample run tags)
rule merge_bams:
    input:
        # Get BAM files from all samples in this run tag group
        bams = lambda wildcards: [get_aligned_bam({"sample": s}) for s in RUN_TAG_SAMPLES.get(wildcards.run_tag, [wildcards.run_tag])],
        # Get alignment flag files for all samples in this run tag group
        alignment_flags = lambda wildcards: [get_alignment_flag(s) for s in RUN_TAG_SAMPLES.get(wildcards.run_tag, [wildcards.run_tag])]
    output:
        bam = get_processing_path("merged/{run_tag}/{run_tag}.bam"),
        bai = get_processing_path("merged/{run_tag}/{run_tag}.bam.bai"),
        flag = get_processing_path("merged/{run_tag}/merge_complete.flag")
    params:
        # Get the number of input BAM files for this run tag
        num_inputs = lambda wildcards: len(RUN_TAG_SAMPLES.get(wildcards.run_tag, [wildcards.run_tag])),
        # Pass the sample IDs for logging
        samples = lambda wildcards: ", ".join(RUN_TAG_SAMPLES.get(wildcards.run_tag, [wildcards.run_tag]))
    log:
        get_processing_path("merged/{run_tag}/logs/merge_bams.log")
    benchmark:
        get_processing_path("merged/{run_tag}/benchmarks/merge_bams.benchmark.txt")
    threads: 
        config["cores_default"]
    resources:
        mem_mb = 8000
    conda:
        "../../envs/samtools.yaml"
    singularity:
        config.get("singularity_image", "")
    shell:
        """
        # Create output directories
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log})
        
        # Log start
        echo "Merging/copying BAM files for run tag {wildcards.run_tag}" > {log}
        echo "Input BAMs: {input.bams}" >> {log}
        echo "Output BAM: {output.bam}" >> {log}
        echo "Number of input BAMs: {params.num_inputs}" >> {log}
        echo "Samples: {params.samples}" >> {log}
        
        # Use different approaches based on number of input files
        if [ {params.num_inputs} -eq 1 ]; then
            # For single-sample run tags, just copy the file
            echo "Single-sample run tag detected, copying file" >> {log}
            cp {input.bams[0]} {output.bam}
            cp {input.bams[0]}.bai {output.bai}
        else
            # For multi-sample run tags, merge the files
            echo "Multi-sample run tag detected, merging files" >> {log}
            samtools merge -@ {threads} -f {output.bam} {input.bams} 2>> {log}
            
            # Create index for the merged BAM file
            echo "Creating index for merged BAM" >> {log}
            samtools index -@ {threads} {output.bam} 2>> {log}
        fi
        
        # Check if output was created successfully
        if [[ ! -s {output.bam} ]]; then
            echo "ERROR: Failed to create merged BAM file" >> {log}
            exit 1
        fi
        
        if [[ ! -s {output.bai} ]]; then
            echo "ERROR: Failed to create BAM index" >> {log}
            exit 1
        fi
        
        # Create flag file to indicate completion
        echo "BAM merging complete for run tag {wildcards.run_tag}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "Files created:" >> {output.flag}
        echo "- BAM file: {output.bam}" >> {output.flag}
        echo "- BAM index: {output.bai}" >> {output.flag}
        echo "Merged samples: {params.samples}" >> {output.flag}
        
        echo "BAM merge complete for {wildcards.run_tag}" >> {log}
        """