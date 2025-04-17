"""
Rules for copying final results to the results directory
Updated to handle run tag groups with a simpler approach
"""
import os
import shutil

# Rule to copy final results to the results directory
rule copy_results:
    input:
        # Required inputs from previous steps - use helper functions to get type-specific flags
        qc_complete = lambda wildcards: get_qc_flag(wildcards.run_tag),
        featurecounts_complete = lambda wildcards: get_featurecounts_flag(wildcards.run_tag),
        coverage_complete = lambda wildcards: get_coverage_flag(wildcards.run_tag),
        benchmark_summary = get_processing_path("{run_tag}/benchmarks_summary.txt")
    output:
        # Flag to indicate copying is complete
        complete = get_results_path("{run_tag}/copy_complete.txt")
    params:
        # Source directories and files
        proc_dir = get_processing_path("{run_tag}"),
        # BAM file
        bam_file = get_picard_bam,
        bam_index = lambda wildcards: f"{get_picard_bam(wildcards)}.bai",
        # Coverage files
        all_bw = get_processing_path("{run_tag}/{run_tag}_all.bw"),
        unique_bw = get_processing_path("{run_tag}/{run_tag}_unique.bw"),
        # Merged benchmarks file
        merged_benchmarks = get_processing_path("{run_tag}/benchmarks/merged_benchmarks.txt"),
        # QC directory
        qc_dir = get_processing_path("{run_tag}/qc"),
        # Benchmark directory
        benchmark_dir = get_processing_path("{run_tag}/benchmarks"),
        # Logs directory
        logs_dir = get_processing_path("{run_tag}/logs"),
        # Original samples in this run tag (for run tags)
        original_samples = lambda wildcards: RUN_TAG_SAMPLES.get(wildcards.run_tag, [wildcards.run_tag]),
        # Target directories
        results_sample_dir = get_results_path("{run_tag}"),
        results_qc_dir = get_results_path("{run_tag}/qc"),
        results_benchmark_dir = get_results_path("{run_tag}/benchmarks"),
        results_logs_dir = get_results_path("{run_tag}/logs"),
        # Config options for copying
        copy_bam = config.get("copy_bam", False),
        copy_fastq = config.get("copy_fastq", False),
        copy_benchmarks = config.get("copy_benchmarks", True),
        copy_logs = config.get("copy_logs", True)
    log:
        get_processing_path("{run_tag}/logs/copy_results.log")
    shell:
        """
        # Create target directories
        mkdir -p {params.results_sample_dir}
        mkdir -p {params.results_qc_dir}
        mkdir -p $(dirname {log})
        
        # Start logging
        echo "Copying results for {wildcards.run_tag} to {params.results_sample_dir}" > {log}
        echo "Timestamp: $(date)" >> {log}
        
        # Log original samples if this is a run tag group
        if [ "{wildcards.run_tag}" != "{params.original_samples[0]}" ]; then
            echo "Run tag group containing samples: {params.original_samples}" >> {log}
        fi
        
        # Copy QC directory
        echo "Copying QC directory" >> {log}
        cp -R {params.qc_dir}/* {params.results_qc_dir}/ 2>> {log}
        
        # Copy bigWig files
        echo "Copying bigWig files" >> {log}
        cp {params.all_bw} {params.results_sample_dir}/ 2>> {log}
        cp {params.unique_bw} {params.results_sample_dir}/ 2>> {log}
        
        # Copy count files - use the naming pattern for both paired and single-end
        echo "Copying count files" >> {log}
        cp {params.proc_dir}/{wildcards.run_tag}_counts_*.txt {params.results_sample_dir}/ 2>> {log}
        
        # Copy benchmark summary file
        echo "Copying benchmark summary" >> {log}
        cp {params.merged_benchmarks} {params.results_sample_dir}/ 2>> {log}
        
        # Copy BAM file if configured - handle both capitalizations
        if [ "{params.copy_bam}" = "True" ] || [ "{params.copy_bam}" = "true" ]; then
            echo "Copying BAM file and index" >> {log}
            cp {params.bam_file} {params.results_sample_dir}/ 2>> {log}
            cp {params.bam_index} {params.results_sample_dir}/ 2>> {log}
        else
            echo "Skipping BAM file copy (disabled in config)" >> {log}
        fi
        
        # Copy FASTQ files if configured - handle both capitalizations
        if [ "{params.copy_fastq}" = "True" ] || [ "{params.copy_fastq}" = "true" ]; then
            if [ "{wildcards.run_tag}" != "{params.original_samples[0]}" ]; then
                # For run tag groups, create a subdirectory with original FASTQ files
                mkdir -p {params.results_sample_dir}/original_fastq
                echo "Copying original FASTQ files to subdirectory" >> {log}
                
                # Copy original FASTQ files for each sample in the run tag group
                for sample in {params.original_samples}; do
                    # Determine if paired or single end
                    if grep -q "paired" {params.proc_dir}/qc_paired_complete.flag 2>/dev/null; then
                        if [ -f $(dirname {params.proc_dir})/$sample/$sample.1.fastq.gz ]; then
                            echo "Copying paired-end FASTQ for $sample" >> {log}
                            cp $(dirname {params.proc_dir})/$sample/$sample.1.fastq.gz {params.results_sample_dir}/original_fastq/ 2>> {log}
                            cp $(dirname {params.proc_dir})/$sample/$sample.2.fastq.gz {params.results_sample_dir}/original_fastq/ 2>> {log}
                        fi
                    else
                        if [ -f $(dirname {params.proc_dir})/$sample/$sample.fastq.gz ]; then
                            echo "Copying single-end FASTQ for $sample" >> {log}
                            cp $(dirname {params.proc_dir})/$sample/$sample.fastq.gz {params.results_sample_dir}/original_fastq/ 2>> {log}
                        fi
                    fi
                done
            else
                # For standalone samples, copy FASTQ files directly
                echo "Copying FASTQ files" >> {log}
                if grep -q "paired" {params.proc_dir}/qc_paired_complete.flag 2>/dev/null; then
                    if [ -f $(dirname {params.proc_dir})/{wildcards.run_tag}/{wildcards.run_tag}.1.fastq.gz ]; then
                        cp $(dirname {params.proc_dir})/{wildcards.run_tag}/{wildcards.run_tag}.1.fastq.gz {params.results_sample_dir}/ 2>> {log}
                        cp $(dirname {params.proc_dir})/{wildcards.run_tag}/{wildcards.run_tag}.2.fastq.gz {params.results_sample_dir}/ 2>> {log}
                    fi
                else
                    if [ -f $(dirname {params.proc_dir})/{wildcards.run_tag}/{wildcards.run_tag}.fastq.gz ]; then
                        cp $(dirname {params.proc_dir})/{wildcards.run_tag}/{wildcards.run_tag}.fastq.gz {params.results_sample_dir}/ 2>> {log}
                    fi
                fi
            fi
        else
            echo "Skipping FASTQ file copy (disabled in config)" >> {log}
        fi
        
        # Copy benchmarks directory if configured - handle both capitalizations
        if [ "{params.copy_benchmarks}" = "True" ] || [ "{params.copy_benchmarks}" = "true" ]; then
            echo "Copying benchmarks directory" >> {log}
            mkdir -p {params.results_benchmark_dir}
            cp -R {params.benchmark_dir}/* {params.results_benchmark_dir}/ 2>> {log}
        else
            echo "Skipping benchmarks directory copy (disabled in config)" >> {log}
        fi
        
        # Copy logs directory if configured - handle both capitalizations
        if [ "{params.copy_logs}" = "True" ] || [ "{params.copy_logs}" = "true" ]; then
            echo "Copying logs directory" >> {log}
            mkdir -p {params.results_logs_dir}
            cp -R {params.logs_dir}/* {params.results_logs_dir}/ 2>> {log}
        else
            echo "Skipping logs directory copy (disabled in config)" >> {log}
        fi
        
        # Create completion file with summary of copied files
        echo "Results copied successfully for {wildcards.run_tag}" > {output.complete}
        echo "Timestamp: $(date)" >> {output.complete}
        echo "----------------------------------------" >> {output.complete}
        
        # Add run tag info if applicable
        if [ "{wildcards.run_tag}" != "{params.original_samples[0]}" ]; then
            echo "Run tag group containing samples: {params.original_samples}" >> {output.complete}
            echo "----------------------------------------" >> {output.complete}
        fi
        
        echo "Copied:" >> {output.complete}
        echo "- QC results" >> {output.complete}
        echo "- Coverage tracks (bigWig)" >> {output.complete}
        echo "- Feature counts" >> {output.complete}
        echo "- Benchmark summary" >> {output.complete}
        
        # Also handle both capitalizations in the report
        if [ "{params.copy_bam}" = "True" ] || [ "{params.copy_bam}" = "true" ]; then
            echo "- BAM file and index" >> {output.complete}
        fi
        
        if [ "{params.copy_fastq}" = "True" ] || [ "{params.copy_fastq}" = "true" ]; then
            echo "- FASTQ files" >> {output.complete}
        fi
        
        if [ "{params.copy_benchmarks}" = "True" ] || [ "{params.copy_benchmarks}" = "true" ]; then
            echo "- Benchmark files" >> {output.complete}
        fi
        
        if [ "{params.copy_logs}" = "True" ] || [ "{params.copy_logs}" = "true" ]; then
            echo "- Log files" >> {output.complete}
        fi
        
        echo "Results copying completed successfully" >> {log}
        """

# Rule to mark completion of copying results for all samples
rule copy_results_all:
    input:
        sample_copies = expand(get_results_path("{run_tag}/copy_complete.txt"), run_tag=EFFECTIVE_SAMPLES_TO_PROCESS),
        benchmark_report = get_processing_path("benchmarks_project_summary.txt")
    output:
        project_complete = get_results_path("copy_complete_all.txt")
        # Removed entities_list output
    params:
        benchmark_dest = get_results_path("benchmarks_project_summary.txt"),
        num_samples = len(SAMPLES),
        num_run_tags = len(EFFECTIVE_SAMPLES_TO_PROCESS),
        copy_fastq = config.get("copy_fastq", False),
        copy_bam = config.get("copy_bam", False)
    log:
        get_results_path("logs/copy_results_all.log")
    run:
        # Python code block instead of shell to avoid complex shell scripting
        import os
        from datetime import datetime
        
        # Create log directory
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        
        # Copy the project-wide benchmark summary
        shutil.copy(input.benchmark_report, params.benchmark_dest)
        
        # Create completion file - ONLY KEEP THE FIRST TWO LINES
        with open(output.project_complete, 'w') as out_file:
            out_file.write(f"All results have been copied to the results directory\n")
            out_file.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            # All other content has been removed
        
        # Write to log
        with open(log[0], 'w') as log_file:
            log_file.write("Completed copying all results\n")

# Rule for cleaning up processing directory after results are copied (optional)
rule cleanup_processing:
    input:
        copy_complete = get_results_path("copy_complete_all.txt")
    output:
        cleanup_complete = get_results_path("cleanup_complete.txt")
    params:
        proc_dir = config['processing_dir'],
        cleanup_processing = config.get("cleanup_processing", False)
    log:
        get_results_path("logs/cleanup.log")
    shell:
        """
        # Create logs directory in results
        mkdir -p $(dirname {log})
        
        # Start logging
        echo "Cleanup processing" > {log}
        echo "Timestamp: $(date)" >> {log}
        
        # Handle both capitalizations
        if [ "{params.cleanup_processing}" = "True" ] || [ "{params.cleanup_processing}" = "true" ]; then
            echo "Removing processing directory: {params.proc_dir}" >> {log}
            rm -rf {params.proc_dir} 2>> {log}
            echo "Processing directory removed" >> {log}
            
            echo "Processing directory has been cleaned up" > {output.cleanup_complete}
            echo "Timestamp: $(date)" >> {output.cleanup_complete}
        else
            echo "Cleanup is disabled in config (cleanup_processing=False), keeping processing directory" >> {log}
            echo "Directory preserved: {params.proc_dir}" >> {log}
            
            echo "Cleanup was skipped (disabled in configuration)" > {output.cleanup_complete}
            echo "Timestamp: $(date)" >> {output.cleanup_complete}
            echo "Processing directory has been preserved: {params.proc_dir}" >> {output.cleanup_complete}
        fi
        
        echo "Cleanup step completed" >> {log}
        """