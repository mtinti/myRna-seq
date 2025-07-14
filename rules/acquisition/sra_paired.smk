"""
Rules for downloading paired-end reads directly from NCBI SRA
"""
import os
import re

# Download FASTQ files from SRA
rule download_sra:
    wildcard_constraints:
        sample = "|".join([re.escape(s) for s in SRA_PAIRED_SAMPLES])
    output:
        r1 = temp(get_processing_path("{sample}/sra_download/{sample}_1.fastq")),
        r2 = temp(get_processing_path("{sample}/sra_download/{sample}_2.fastq"))
    params:
        accession = lambda w: SAMPLES_DF.loc[w.sample, 'file_path_1'],
        outdir = lambda w: get_processing_path(f"{w.sample}/sra_download"),
        tmpdir = lambda w: get_processing_path(f"{w.sample}/sra_tmp")
    log:
        get_processing_path("{sample}/logs/sra_download.log")
    threads: 8
    retries: 3
    conda:
        "../../envs/sra.yaml"
    shell:
        """
        mkdir -p {params.outdir} {params.tmpdir}
        fasterq-dump --threads {threads} --split-files \
            --outdir {params.outdir} --temp {params.tmpdir} {params.accession}
        mv {params.outdir}/{params.accession}_1.fastq {output.r1}
        mv {params.outdir}/{params.accession}_2.fastq {output.r2}
        """

# Validate SRA download
rule validate_sra_download:
    input:
        r1 = get_processing_path("{sample}/sra_download/{sample}_1.fastq"),
        r2 = get_processing_path("{sample}/sra_download/{sample}_2.fastq")
    output:
        flag = get_processing_path("{sample}/sra_download/.validated")
    run:
        if os.path.getsize(input.r1) == 0 or os.path.getsize(input.r2) == 0:
            raise ValueError(f"Empty files downloaded for {wildcards.sample}")
        with open(input.r1) as f1, open(input.r2) as f2:
            lines1 = sum(1 for _ in f1)
            lines2 = sum(1 for _ in f2)
            if lines1 != lines2:
                raise ValueError(f"Unequal read counts: R1={lines1/4}, R2={lines2/4}")
        open(output.flag, 'w').close()

# Reformat headers
rule reformat_sra_headers:
    input:
        r1 = get_processing_path("{sample}/sra_download/{sample}_1.fastq"),
        r2 = get_processing_path("{sample}/sra_download/{sample}_2.fastq"),
        valid = get_processing_path("{sample}/sra_download/.validated")
    output:
        r1 = temp(get_processing_path("{sample}/sra_reformatted/{sample}_1.fq")),
        r2 = temp(get_processing_path("{sample}/sra_reformatted/{sample}_2.fq"))
    shell:
        """
        awk '{{print (NR%4==1) ? "@1_" ++i " READ/1" : $0}}' {input.r1} > {output.r1}
        awk '{{print (NR%4==1) ? "@1_" ++i " READ/2" : $0}}' {input.r2} > {output.r2}
        """

# Compress reformatted files
rule gzip_sra_fastq:
    input:
        r1 = get_processing_path("{sample}/sra_reformatted/{sample}_1.fq"),
        r2 = get_processing_path("{sample}/sra_reformatted/{sample}_2.fq")
    output:
        r1 = get_processing_path("{sample}/{sample}_R1.fastq.gz"),
        r2 = get_processing_path("{sample}/{sample}_R2.fastq.gz")
    threads: 4
    conda:
        "../../envs/core.yaml"
    shell:
        """
        gzip -c {input.r1} > {output.r1}
        gzip -c {input.r2} > {output.r2}
        """
