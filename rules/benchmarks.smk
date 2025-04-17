"""
Rules for merging all benchmark files into a single report
Updated to work with run tags - with Python-based merge_benchmarks_all
"""
import os
import glob
import re

# Rule to merge all benchmark files for a run tag
rule merge_benchmarks:
    input:
        # Depend on the completion of previous steps
        qc_complete = lambda wildcards: get_qc_flag(wildcards.run_tag),
        coverage_complete = lambda wildcards: get_coverage_flag(wildcards.run_tag),
        featurecounts_complete = lambda wildcards: get_featurecounts_flag(wildcards.run_tag)
    output:
        merged = get_processing_path("{run_tag}/benchmarks/merged_benchmarks.txt")
    log:
        get_processing_path("{run_tag}/logs/merge_benchmarks.log")
    params:
        benchmark_dir = get_processing_path("{run_tag}/benchmarks")
    resources:
        mem_mb = lambda wildcards: config.get("mem_merge", 1000)
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.merged})
        
        # Start logging
        echo "Merging benchmark files for {wildcards.run_tag}" > {log}
        echo "Output: {output.merged}" >> {log}
        
        # Write header with added 'step' column
        echo -e "step\\ts\\th:m:s\\tmax_rss\\tmax_vms\\tmax_uss\\tmax_pss\\tio_in\\tio_out\\tmean_load\\tcpu_time" > {output.merged}
        
        # Process each benchmark file in the benchmarks directory
        find {params.benchmark_dir} -name "*.benchmark.txt" | sort | while read file; do
            # Extract the step name (text before first dot in filename)
            fname=$(basename "$file")
            step=$(echo "$fname" | cut -d. -f1)
            
            # Skip the first line (header) and add step to each line
            tail -n +2 "$file" | sed "s/^/$step\\t/" >> {output.merged} 2>> {log}
            
            echo "Processed: $file -> step: $step" >> {log}
        done
        
        # Count processed files
        count=$(grep -c "^" {output.merged})
        count=$((count - 1))  # Subtract 1 for the header line
        echo "Merged $count benchmark entries successfully" >> {log}
        
        # Check if merged file contains data beyond the header
        if [ "$count" -lt 1 ]; then
            echo "WARNING: No benchmark data was found to merge" >> {log}
            # Don't fail the pipeline, just create a note in the file
            echo "No benchmark data was found" >> {output.merged}
        else
            echo "Benchmark files merged successfully" >> {log}
        fi
        """

# Rule for creating a run tag-level summary of all benchmarks
rule merge_benchmarks_summary:
    input:
        merged = get_processing_path("{run_tag}/benchmarks/merged_benchmarks.txt")
    output:
        summary = get_processing_path("{run_tag}/benchmarks_summary.txt")
    log:
        get_processing_path("{run_tag}/logs/benchmarks_summary.log")
    shell:
        """
        # Create log directory
        mkdir -p $(dirname {log})
        
        # Start logging
        echo "Creating benchmark summary for {wildcards.run_tag}" > {log}
        
        # Extract key statistics from the merged benchmarks
        echo "BENCHMARK SUMMARY FOR {wildcards.run_tag}" > {output.summary}
        echo "Generated: $(date)" >> {output.summary}
        echo "----------------------------------------" >> {output.summary}
        
        # Count steps
        steps=$(tail -n +2 {input.merged} | cut -f1 | sort | uniq | wc -l)
        echo "Number of steps: $steps" >> {output.summary}
        
        # Total runtime
        total_s=$(tail -n +2 {input.merged} | awk '{{sum+=$2}} END {{printf "%.2f", sum}}')
        echo "Total runtime (seconds): $total_s" >> {output.summary}
        
        # Maximum memory
        max_mem=$(tail -n +2 {input.merged} | sort -k4,4nr | head -n 1)
        if [ -n "$max_mem" ]; then
            step=$(echo "$max_mem" | cut -f1)
            mem=$(echo "$max_mem" | cut -f4)
            echo "Maximum memory (RSS): $mem (step: $step)" >> {output.summary}
        else
            echo "Maximum memory (RSS): No data" >> {output.summary}
        fi
        
        # Maximum CPU time
        max_cpu=$(tail -n +2 {input.merged} | sort -k10,10nr | head -n 1)
        if [ -n "$max_cpu" ]; then
            step=$(echo "$max_cpu" | cut -f1)
            cpu=$(echo "$max_cpu" | cut -f10)
            echo "Maximum CPU time: $cpu (step: $step)" >> {output.summary}
        else
            echo "Maximum CPU time: No data" >> {output.summary}
        fi
        
        echo "----------------------------------------" >> {output.summary}
        echo "See {input.merged} for full details" >> {output.summary}
        
        echo "Benchmark summary created successfully" >> {log}
        """

# Rule for combining all run tag benchmark summaries into a project-wide report - MODIFIED
rule merge_benchmarks_all:
    input:
        summaries = expand(get_processing_path("{run_tag}/benchmarks_summary.txt"), run_tag=EFFECTIVE_SAMPLES_TO_PROCESS)
    output:
        project_report = get_processing_path("benchmarks_project_summary.txt")
    params:
        # Pre-calculate the number of samples and run tags
        num_samples = len(SAMPLES),
        num_run_tags = len(EFFECTIVE_SAMPLES_TO_PROCESS)
    log:
        get_processing_path("logs/benchmarks_project_summary.log")
    run:
        # Python code block instead of shell to avoid complex shell scripting
        import os
        from datetime import datetime
        
        # Create log directory
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        
        # Start logging
        with open(log[0], 'w') as log_file:
            log_file.write("Creating project-wide benchmark summary\n")
        
        # Create header for project report - ONLY KEEP THE FIRST TWO LINES
        with open(output.project_report, 'w') as report:
            report.write("PROJECT-WIDE BENCHMARK SUMMARY\n")
            report.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            # All other content has been removed
        
        # Append to log that we finished
        with open(log[0], 'a') as log_file:
            log_file.write("Project benchmark summary created successfully\n")