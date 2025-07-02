"""
Rules for copying final results to the results directory
Updated to handle run tag groups with a simpler approach
"""
import os
import shutil

# Rule to copy final results to the results directory
# Rule to copy final results to the results directory and optionally clean up processing files
# Rule to copy final results to the results directory and optionally clean up processing files
# Rule to copy final results to the results directory and optionally clean up processing files
rule copy_results:
    input:
        # Required inputs from previous steps - use helper functions to get type-specific flags
        qc_complete = lambda wildcards: get_qc_flag(wildcards.run_tag),
        featurecounts_complete = lambda wildcards: get_featurecounts_flag(wildcards.run_tag),
        coverage_complete = lambda wildcards: get_coverage_flag(wildcards.run_tag),
        benchmark_summary = get_processing_path("{run_tag}/benchmarks_summary.txt")
    output:
        # Flag in results directory (won't be cleaned up)
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
        
        # Config options for copying - convert to lowercase for consistency
        copy_bam = str(config.get("copy_bam", False)).lower(),
        copy_fastq = str(config.get("copy_fastq", False)).lower(),
        copy_benchmarks = str(config.get("copy_benchmarks", True)).lower(),
        copy_logs = str(config.get("copy_logs", True)).lower(),
        
        # Cleanup option
        cleanup_processing = str(config.get("cleanup_processing", False)).lower(),
        
        # Directories to clean up after copying
        sample_dirs_to_cleanup = lambda wildcards: [
            get_processing_path(sample) 
            for sample in RUN_TAG_SAMPLES.get(wildcards.run_tag, [wildcards.run_tag])
        ],
        merged_dir_to_cleanup = get_processing_path("merged/{run_tag}"),
        run_tag_dir_to_cleanup = get_processing_path("{run_tag}")
        
    log:
        # Log file in results directory (won't be cleaned up)
        get_results_path("{run_tag}/logs/copy_results.log")
    shell:
        """
        # Create target directories
        mkdir -p {params.results_sample_dir}
        mkdir -p {params.results_qc_dir}
        mkdir -p $(dirname {log})
        
        # Start logging
        echo "Copying results for {wildcards.run_tag} to {params.results_sample_dir}" > {log}
        echo "Timestamp: $(date)" >> {log}
        echo "Cleanup setting: {params.cleanup_processing}" >> {log}
        
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
        
        # Copy count files
        echo "Copying count files" >> {log}
        cp {params.proc_dir}/{wildcards.run_tag}_counts_*.txt {params.results_sample_dir}/ 2>> {log}
        
        # Copy benchmark summary file
        echo "Copying benchmark summary" >> {log}
        cp {params.merged_benchmarks} {params.results_sample_dir}/ 2>> {log}
        
        # Copy BAM file if configured
        if [ "{params.copy_bam}" = "true" ]; then
            echo "Copying BAM file and index" >> {log}
            cp {params.bam_file} {params.results_sample_dir}/ 2>> {log}
            cp {params.bam_index} {params.results_sample_dir}/ 2>> {log}
        else
            echo "Skipping BAM file copy (disabled in config)" >> {log}
        fi
        
        # Copy FASTQ files if configured
        if [ "{params.copy_fastq}" = "true" ]; then
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
        
        # Copy benchmarks directory if configured
        if [ "{params.copy_benchmarks}" = "true" ]; then
            echo "Copying benchmarks directory" >> {log}
            mkdir -p {params.results_benchmark_dir}
            cp -R {params.benchmark_dir}/* {params.results_benchmark_dir}/ 2>> {log}
        else
            echo "Skipping benchmarks directory copy (disabled in config)" >> {log}
        fi
        
        # Copy logs directory if configured
        if [ "{params.copy_logs}" = "true" ]; then
            echo "Copying logs directory" >> {log}
            mkdir -p {params.results_logs_dir}
            cp -R {params.logs_dir}/* {params.results_logs_dir}/ 2>> {log}
        else
            echo "Skipping logs directory copy (disabled in config)" >> {log}
        fi
        
        # Create completion file BEFORE cleanup (so it exists even if cleanup fails)
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
        
        if [ "{params.copy_bam}" = "true" ]; then
            echo "- BAM file and index" >> {output.complete}
        fi
        
        if [ "{params.copy_fastq}" = "true" ]; then
            echo "- FASTQ files" >> {output.complete}
        fi
        
        if [ "{params.copy_benchmarks}" = "true" ]; then
            echo "- Benchmark files" >> {output.complete}
        fi
        
        if [ "{params.copy_logs}" = "true" ]; then
            echo "- Log files" >> {output.complete}
        fi
        
        # Clean up processing directories if configured (AFTER creating completion file)
        if [ "{params.cleanup_processing}" = "true" ]; then
            echo "----------------------------------------" >> {log}
            echo "Starting cleanup of processing directories..." >> {log}
            
            # Wait a moment to ensure all file operations are complete
            sleep 1
            
            # Clean up individual sample directories
            for sample_dir in {params.sample_dirs_to_cleanup}; do
                if [ -d "$sample_dir" ]; then
                    echo "Removing sample directory: $sample_dir" >> {log}
                    
                    # First remove all contents recursively
                    rm -rf "$sample_dir"/* 2>> {log}
                    rm -rf "$sample_dir"/.[!.]* 2>> {log}  # Remove hidden files/dirs
                    
                    # Then remove the directory itself
                    rmdir "$sample_dir" 2>> {log}
                    
                    # Force removal if rmdir failed
                    if [ -d "$sample_dir" ]; then
                        echo "Directory still exists, forcing removal..." >> {log}
                        rm -rf "$sample_dir" 2>> {log}
                    fi
                    
                    if [ ! -d "$sample_dir" ]; then
                        echo "Successfully removed: $sample_dir" >> {log}
                    else
                        echo "WARNING: Failed to completely remove: $sample_dir" >> {log}
                        ls -la "$sample_dir" 2>> {log} || echo "Cannot list directory contents" >> {log}
                    fi
                else
                    echo "Sample directory does not exist: $sample_dir" >> {log}
                fi
            done
            
            # Clean up merged directory for this run tag
            if [ -d "{params.merged_dir_to_cleanup}" ]; then
                echo "Removing merged directory: {params.merged_dir_to_cleanup}" >> {log}
                
                # Remove contents first
                rm -rf "{params.merged_dir_to_cleanup}"/* 2>> {log}
                rm -rf "{params.merged_dir_to_cleanup}"/.[!.]* 2>> {log}
                
                # Then remove the directory
                rmdir "{params.merged_dir_to_cleanup}" 2>> {log}
                
                # Force removal if needed
                if [ -d "{params.merged_dir_to_cleanup}" ]; then
                    echo "Directory still exists, forcing removal..." >> {log}
                    rm -rf "{params.merged_dir_to_cleanup}" 2>> {log}
                fi
                
                if [ ! -d "{params.merged_dir_to_cleanup}" ]; then
                    echo "Successfully removed: {params.merged_dir_to_cleanup}" >> {log}
                else
                    echo "WARNING: Failed to completely remove: {params.merged_dir_to_cleanup}" >> {log}
                    ls -la "{params.merged_dir_to_cleanup}" 2>> {log} || echo "Cannot list directory contents" >> {log}
                fi
            else
                echo "Merged directory does not exist: {params.merged_dir_to_cleanup}" >> {log}
            fi
            
            # Clean up run tag directory
            if [ -d "{params.run_tag_dir_to_cleanup}" ]; then
                echo "Removing run tag directory: {params.run_tag_dir_to_cleanup}" >> {log}
                
                # Remove contents first
                rm -rf "{params.run_tag_dir_to_cleanup}"/* 2>> {log}
                rm -rf "{params.run_tag_dir_to_cleanup}"/.[!.]* 2>> {log}
                
                # Then remove the directory
                rmdir "{params.run_tag_dir_to_cleanup}" 2>> {log}
                
                # Force removal if needed
                if [ -d "{params.run_tag_dir_to_cleanup}" ]; then
                    echo "Directory still exists, forcing removal..." >> {log}
                    rm -rf "{params.run_tag_dir_to_cleanup}" 2>> {log}
                fi
                
                # Also try to remove the parent merged directory if it's empty
                merged_parent=$(dirname "{params.merged_dir_to_cleanup}")
                if [ -d "$merged_parent" ] && [ -z "$(ls -A "$merged_parent" 2>/dev/null)" ]; then
                    echo "Removing empty merged parent directory: $merged_parent" >> {log}
                    rmdir "$merged_parent" 2>> {log}
                fi
                
                if [ ! -d "{params.run_tag_dir_to_cleanup}" ]; then
                    echo "Successfully removed: {params.run_tag_dir_to_cleanup}" >> {log}
                else
                    echo "WARNING: Failed to completely remove: {params.run_tag_dir_to_cleanup}" >> {log}
                    ls -la "{params.run_tag_dir_to_cleanup}" 2>> {log} || echo "Cannot list directory contents" >> {log}
                fi
            else
                echo "Run tag directory does not exist: {params.run_tag_dir_to_cleanup}" >> {log}
            fi
            
            # Final verification and reporting
            echo "Cleanup verification:" >> {log}
            remaining_dirs=""
            for sample_dir in {params.sample_dirs_to_cleanup}; do
                if [ -d "$sample_dir" ]; then
                    remaining_dirs="$remaining_dirs $sample_dir"
                fi
            done
            if [ -d "{params.merged_dir_to_cleanup}" ]; then
                remaining_dirs="$remaining_dirs {params.merged_dir_to_cleanup}"
            fi
            if [ -d "{params.run_tag_dir_to_cleanup}" ]; then
                remaining_dirs="$remaining_dirs {params.run_tag_dir_to_cleanup}"
            fi
            
            if [ -z "$remaining_dirs" ]; then
                echo "All directories successfully removed" >> {log}
            else
                echo "WARNING: Some directories could not be removed:$remaining_dirs" >> {log}
            fi
            
            echo "Cleanup completed for {wildcards.run_tag}" >> {log}
            echo "----------------------------------------" >> {log}
            
            # Add cleanup info to completion file
            echo "----------------------------------------" >> {output.complete}
            echo "CLEANUP PERFORMED:" >> {output.complete}
            echo "Removed processing directories:" >> {output.complete}
            for sample_dir in {params.sample_dirs_to_cleanup}; do
                echo "- $sample_dir" >> {output.complete}
            done
            echo "- {params.merged_dir_to_cleanup}" >> {output.complete}
            echo "- {params.run_tag_dir_to_cleanup}" >> {output.complete}
        else
            echo "Cleanup is disabled, preserving processing directories" >> {log}
        fi
        
        echo "Results copying and cleanup completed successfully" >> {log}
        """
        
# Rule to mark completion of copying results for all samples
rule copy_results_all:
    input:
        sample_copies = expand(get_results_path("{run_tag}/copy_complete.txt"), run_tag=EFFECTIVE_SAMPLES_TO_PROCESS),
        benchmark_report = get_processing_path("benchmarks_project_summary.txt")
    output:
        project_complete = get_results_path("copy_complete_all.txt")
    params:
        benchmark_dest = get_results_path("benchmarks_project_summary.txt"),
        num_samples = len(SAMPLES),
        num_run_tags = len(EFFECTIVE_SAMPLES_TO_PROCESS),
        copy_fastq = config.get("copy_fastq", False),
        copy_bam = config.get("copy_bam", False),
        cleanup_processing = config.get("cleanup_processing", False)
    log:
        get_results_path("logs/copy_results_all.log")
    shell:
        """
        # Create log directory
        mkdir -p $(dirname {log})
        
        # Copy the project-wide benchmark summary
        cp {input.benchmark_report} {params.benchmark_dest} 2>> {log}
        
        # Create completion file
        echo "All results have been copied to the results directory" > {output.project_complete}
        echo "Timestamp: $(date)" >> {output.project_complete}
        echo "Number of samples processed: {params.num_samples}" >> {output.project_complete}
        echo "Number of run tags/groups: {params.num_run_tags}" >> {output.project_complete}
        echo "Configuration used:" >> {output.project_complete}
        echo "- Copy FASTQ files: {params.copy_fastq}" >> {output.project_complete}
        echo "- Copy BAM files: {params.copy_bam}" >> {output.project_complete}
        echo "- Cleanup processing: {params.cleanup_processing}" >> {output.project_complete}
        
        # Write to log
        echo "Completed copying all results" > {log}
        echo "Timestamp: $(date)" >> {log}
        echo "All individual copy operations completed successfully" >> {log}
        """