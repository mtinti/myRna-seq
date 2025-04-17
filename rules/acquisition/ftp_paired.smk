"""
Rules for handling FTP paired-end input files
"""
import os
import re

rule download_ftp_paired_files:
    wildcard_constraints:
        sample = "|".join([re.escape(s) for s in PAIRED_FTP_SAMPLES])
    output:
        r1 = get_processing_path("{sample}/{sample}.1.fastq.gz"),
        r2 = get_processing_path("{sample}/{sample}.2.fastq.gz"),
        flag = get_processing_path("{sample}/acquisition_paired_complete.flag")
    params:
        url_r1 = lambda wildcards: SAMPLES_DF.loc[wildcards.sample, 'file_path_1'],
        url_r2 = lambda wildcards: SAMPLES_DF.loc[wildcards.sample, 'file_path_2']
    log:
        get_processing_path("{sample}/logs/download_paired_ftp.log")
    benchmark:
        get_processing_path("{sample}/benchmarks/download_paired_ftp.benchmark.txt")
    resources:
        network=1,
        mem_mb = lambda wildcards: config.get("mem_download", 2000)
    conda:
        "../../envs/wget.yaml"
    singularity:
        config.get("singularity_image", "")
    shell:
        """
        # Create output directories
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.r1})
        
        echo "Downloading FTP paired-end files for {wildcards.sample}" > {log}
        
        # Check if files already exist
        if [[ -s {output.r1} && -s {output.r2} ]]; then
            echo "Files already exist, skipping download" >> {log}
        else
            echo "Downloading {params.url_r1} to {output.r1}" >> {log}
            wget -t 5 -nv -c -T 60 -O {output.r1} {params.url_r1} 2>> {log}
            echo "Downloading {params.url_r2} to {output.r2}" >> {log}
            wget -t 5 -nv -c -T 60 -O {output.r2} {params.url_r2} 2>> {log}
        fi
        
        # Create flag file to indicate completion
        echo "Acquisition complete for paired-end sample {wildcards.sample}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "Files downloaded:" >> {output.flag}
        echo "- {output.r1} from {params.url_r1}" >> {output.flag}
        echo "- {output.r2} from {params.url_r2}" >> {output.flag}
        """