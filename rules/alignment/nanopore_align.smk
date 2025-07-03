"""
Rules for alignment of nanopore reads using minimap2
Skip fastp and mark duplicates.
"""
import os
import re

def get_nanopore_align_input(wildcards):
    return {
        'r': get_nanopore_fastq(wildcards),
        'checksums': get_processing_path(f"{wildcards.sample}/checksums_single_verified.flag")
    }

rule align_nanopore:
    wildcard_constraints:
        sample = "|".join([re.escape(s) for s in NANOPORE_LOCAL_SAMPLES + NANOPORE_FTP_SAMPLES]) if NANOPORE_LOCAL_SAMPLES + NANOPORE_FTP_SAMPLES else "^$"
    input:
        unpack(get_nanopore_align_input)
    output:
        bam = get_processing_path("{sample}/{sample}.bam"),
        bai = get_processing_path("{sample}/{sample}.bam.bai"),
        stats = get_processing_path("{sample}/qc/minimap2/{sample}.minimap2_stats.txt"),
        flag = get_processing_path("{sample}/alignment_single_complete.flag")
    params:
        preset = config.get("nanopore_preset", "map-ont"),
        genome_index = config["processing_genome_index"],
        rg_id = lambda wc: wc.sample,
        rg_sm = lambda wc: wc.sample,
        rg_lb = lambda wc: f"{wc.sample}_lib",
        rg_pu = lambda wc: f"{wc.sample}_unit",
        rg_pl = "ONT"
    log:
        get_processing_path("{sample}/logs/alignment_nanopore.log")
    benchmark:
        get_processing_path("{sample}/benchmarks/align_nanopore.benchmark.txt")
    threads:
        config.get("cores_nanopore_align", 8)
    resources:
        mem_mb = config.get("mem_nanopore_align", 16000)
    conda:
        "../../envs/nanopore_alignment.yaml"
    singularity:
        config.get("singularity_image", "")
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {output.stats})

        echo "Aligning nanopore reads for {wildcards.sample}" > {log}
        echo "Input FASTQ: {input.r}" >> {log}
        echo "Using preset: {params.preset}" >> {log}
        echo "Genome index base: {params.genome_index}" >> {log}

        minimap2 -ax {params.preset} -t {threads} \
            -R '@RG\tID:{params.rg_id}\tSM:{params.rg_sm}\tLB:{params.rg_lb}\tPU:{params.rg_pu}\tPL:{params.rg_pl}' \
            {params.genome_index} {input.r} 2> {output.stats} | \
            samtools view -bSu -@ {threads} | \
            samtools sort -@ {threads} -o {output.bam} >> {log} 2>&1

        samtools index -@ {threads} {output.bam}

        echo "Alignment complete for {wildcards.sample}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "BAM: {output.bam}" >> {output.flag}
        echo "Stats: {output.stats}" >> {output.flag}
        """
