"""
Rule for building Bowtie2 indexes from the staged reference FASTA.
"""

import os

rule build_bowtie2_index:
    input:
        fasta = get_reference_fasta
    output:
        get_bowtie2_index_files()
    params:
        index_prefix = config["processing_genome_index"],
        output_dir = os.path.dirname(get_bowtie2_index_files()[0]),
        log_dir = os.path.dirname(get_processing_path("reference/logs/build_bowtie2_index.log")),
        benchmark_dir = os.path.dirname(get_processing_path("reference/benchmarks/build_bowtie2_index.benchmark.txt")),
    log:
        get_processing_path("reference/logs/build_bowtie2_index.log")
    benchmark:
        get_processing_path("reference/benchmarks/build_bowtie2_index.benchmark.txt")
    threads:
        config.get("cores_index", config.get("cores_align", 8))
    resources:
        mem_mb = 8000
    conda:
        "../../envs/alignment.yaml"
    singularity:
        config.get("singularity_image", "")
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p {params.benchmark_dir}
        mkdir -p {params.output_dir}

        echo "Building Bowtie2 index from {input.fasta}" > {log}
        echo "Index prefix: {params.index_prefix}" >> {log}
        echo "Using {threads} threads" >> {log}

        if ! command -v bowtie2-build &> /dev/null; then
            echo "ERROR: bowtie2-build command not found in environment" >> {log}
            exit 1
        fi

        missing_files=0
        for idx_file in {output}; do
            if [[ ! -s "$idx_file" ]]; then
                missing_files=1
                break
            fi
        done

        if [[ $missing_files -eq 0 ]]; then
            echo "Existing Bowtie2 index detected, skipping rebuild." >> {log}
        else
            echo "Running bowtie2-build" >> {log}
            bowtie2-build --threads {threads} {input.fasta} {params.index_prefix} >> {log} 2>&1
        fi
        """
