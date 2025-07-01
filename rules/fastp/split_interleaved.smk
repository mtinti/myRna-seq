rule split_interleaved:
    input:
        fastq    = get_fastq_single,
        checksums = get_processing_path("{sample}/checksums_single_verified.flag")
    output:
        # Outputs use the standard *.fastq.gz extension so downstream rules can
        # treat the files just like regular paired-end reads.
        r1       = get_processing_path("{sample}/{sample}.1.fastq.gz"),
        r2       = get_processing_path("{sample}/{sample}.2.fastq.gz"),
        flag     = get_processing_path("{sample}/split_interleaved_complete.flag")
    log:
        get_processing_path("{sample}/logs/split_interleaved.log")
    params:
        script   = "split_interleaved.sh",
        prefix   = get_processing_path("{sample}/{sample}")
    shell:
        """
        # Call the provided splitting script and capture all output to the Snakemake log
        bash {params.script} {input.fastq} {params.prefix} > {log} 2>&1

        # Verify that both files were created and are non‐empty
        if [[ ! -s {output.r1} || ! -s {output.r2} ]]; then
            echo "ERROR: Failed to create output files" >> {log}
            exit 1
        fi

        # Touch the completion flag with metadata
        echo "Interleaved FASTQ splitting completed for {wildcards.sample}" > {output.flag}
        echo "Timestamp: $(date)"                            >> {output.flag}
        echo "Input file: {input.fastq}"                      >> {output.flag}
        echo "Output R1: {output.r1}"                         >> {output.flag}
        echo "Output R2: {output.r2}"                         >> {output.flag}
        """
