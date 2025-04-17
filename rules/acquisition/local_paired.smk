"""
Rules for handling local paired-end input files
"""
import os
import re

# Rule specifically for paired-end samples with local source
rule copy_local_paired_files:
    wildcard_constraints:
        sample = "|".join([re.escape(s) for s in PAIRED_LOCAL_SAMPLES])
    output:
        r1 = get_processing_path("{sample}/{sample}.1.fastq.gz"),
        r2 = get_processing_path("{sample}/{sample}.2.fastq.gz"),
        flag = get_processing_path("{sample}/acquisition_paired_complete.flag")
    params:
        src_r1 = lambda wildcards: SAMPLES_DF.loc[wildcards.sample, 'file_path_1'],
        src_r2 = lambda wildcards: SAMPLES_DF.loc[wildcards.sample, 'file_path_2']
    log:
        get_processing_path("{sample}/logs/copy_paired_local.log")
    benchmark:
        get_processing_path("{sample}/benchmarks/copy_paired_local.benchmark.txt")
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
        mkdir -p $(dirname {output.r1})
        
        echo "Copying local paired-end files for {wildcards.sample}" > {log}
        
        # Check if files already exist
        if [[ -s {output.r1} && -s {output.r2} ]]; then
            echo "Files already exist, skipping copy" >> {log}
        else
            echo "Copying {params.src_r1} to {output.r1}" >> {log}
            cp {params.src_r1} {output.r1}
            echo "Copying {params.src_r2} to {output.r2}" >> {log}
            cp {params.src_r2} {output.r2}
        fi
        
        # Create flag file to indicate completion
        echo "Acquisition complete for paired-end sample {wildcards.sample}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "Files copied:" >> {output.flag}
        echo "- {output.r1}" >> {output.flag}
        echo "- {output.r2}" >> {output.flag}
        """