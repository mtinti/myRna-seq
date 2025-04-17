"""
Rules for handling local single-end input files
"""
import os

# Rule specifically for single-end samples with local source
import re

rule copy_local_single_file:
    wildcard_constraints:
        sample = "|".join([re.escape(s) for s in SINGLE_LOCAL_SAMPLES])
    output:
        r = get_processing_path("{sample}/{sample}.fastq.gz"),
        flag = get_processing_path("{sample}/acquisition_single_complete.flag")
    params:
        # Using file_path_1 for single-end data
        src = lambda wildcards: SAMPLES_DF.loc[wildcards.sample, 'file_path_1']
    log:
        get_processing_path("{sample}/logs/copy_single_local.log")
    benchmark:
        get_processing_path("{sample}/benchmarks/copy_single_local.benchmark.txt")
    resources:
        io=1,
        mem_mb = lambda wildcards: config.get("mem_copy", 2000)
    conda:
        "../../envs/core.yaml"
    singularity:
        config.get("singularity_image", "")
    shell:
        """
        # Create output directories
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.r})
        
        echo "Copying local single-end file for {wildcards.sample}" > {log}
        
        # Check if file already exists
        if [[ -s {output.r} ]]; then
            echo "File already exists, skipping copy" >> {log}
        else
            echo "Copying {params.src} to {output.r}" >> {log}
            cp {params.src} {output.r}
        fi
        
        # Create flag file to indicate completion
        echo "Acquisition complete for single-end sample {wildcards.sample}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "File copied:" >> {output.flag}
        echo "- {output.r}" >> {output.flag}
        """