"""
Rules for RNA-seq QC analysis for paired-end data
Uses qualimap rnaseq with paired-end specific options
"""
import os

# Rule for running qualimap rnaseq on paired-end data
rule qualimap_rnaseq_paired:
    input:
        bam = get_picard_bam,
        bai = lambda wildcards: f"{get_picard_bam(wildcards)}.bai",
        bamqc_flag = get_processing_path("{sample}/bamqc_complete.flag")
    output:
        dir = directory(get_processing_path("{sample}/qc/{sample}_qualimap_rnaseq/")),
        flag = get_processing_path("{sample}/qc_paired_complete.flag")
    params:
        memory = config.get("qualimap_memory", "10G"),
        gtf = config["processing_gtf_file"]
    log:
        stderr = get_processing_path("{sample}/logs/qualimap_rnaseq.log"),
        stdout = get_processing_path("{sample}/logs/qualimap_rnaseq.stdout.log")
    benchmark:
        get_processing_path("{sample}/benchmarks/qualimap_rnaseq.benchmark.txt")
    threads: 
        config["cores_default"]
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
        echo "Running qualimap rnaseq (paired-end) on {wildcards.sample}" > {log.stderr}
        echo "Input BAM: {input.bam}" >> {log.stderr}
        echo "GTF file: {params.gtf}" >> {log.stderr}
        echo "Output directory: {output.dir}" >> {log.stderr}
        echo "Java memory: {params.memory}" >> {log.stderr}
        
        # Check if GTF file exists
        if [[ ! -f {params.gtf} ]]; then
            echo "ERROR: GTF file {params.gtf} not found" >> {log.stderr}
            exit 1
        fi
        
        # Run qualimap rnaseq with paired-end option
        qualimap rnaseq --java-mem-size={params.memory} \\
            -bam {input.bam} \\
            -gtf {params.gtf} \\
            -outdir {output.dir} \\
            -outformat HTML \\
            -pe \\
            > {log.stdout} 2>> {log.stderr}
        
        # Check for successful completion by looking for report.html
        if [[ ! -f {output.dir}/qualimapReport.html ]]; then
            echo "ERROR: qualimap rnaseq did not produce report.html" >> {log.stderr}
            exit 1
        fi
        
        # Create flag file to indicate completion
        echo "QC processing completed for paired-end sample {wildcards.sample}" > {output.flag}
        echo "Timestamp: $(date)" >> {output.flag}
        echo "---------------------------------" >> {output.flag}
        echo "RNA-seq QC report: {output.dir}/qualimapReport.html" >> {output.flag}
        
        echo "Qualimap rnaseq (paired-end) completed successfully" >> {log.stderr}
        """