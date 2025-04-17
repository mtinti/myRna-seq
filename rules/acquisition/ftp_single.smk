"""
Rules for handling FTP single-end input files
"""
import os

import re

rule download_ftp_single_file:
    wildcard_constraints:
        sample = "|".join([re.escape(s) for s in SINGLE_FTP_SAMPLES])
    output:
        r = get_processing_path("{sample}/{sample}.fastq.gz"),
        flag = get_processing_path("{sample}/acquisition_single_complete.flag")
    params:
        # Using file_path_1 for single-end data
        url = lambda wildcards: SAMPLES_DF.loc[wildcards.sample, 'file_path_1']
    log:
        get_processing_path("{sample}/logs/download_single_ftp.log")
    benchmark:
        get_processing_path("{sample}/benchmarks/download_single_ftp.benchmark.txt")
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
        mkdir -p $(dirname {output.r})
        
        echo "Downloading FTP single-end file for {wildcards.sample}" > {log}
        
        # Check if file already exists
        if [[ -s {output.r} ]]; then
            echo "File already exists, skipping download" >> {log}
        else
            echo "Downloading {params.url} to {output.r}" >> {log}
            wget -t 5 -nv -c -T 60 -O {output.r} {params.url} 2>> {log}
        fi
        
        # Create flag file to indicate completion
        echo "Acquisition complete for single-end sample {wildcards.sample}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "File downloaded:" >> {output.flag}
        echo "- {output.r} from {params.url}" >> {output.flag}
        """