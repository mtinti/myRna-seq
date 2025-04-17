"""
Utility functions for sample handling in the RNA-seq pipeline.
Extracts common functionality from Snakefile for better organization.
"""
import os
import pandas as pd
import shutil

def ensure_directories_exist(config):
    """Create the necessary directory structure for the pipeline"""
    # Processing directory
    os.makedirs(config['processing_dir'], exist_ok=True)
    
    # Results directory
    os.makedirs(config['results_dir'], exist_ok=True)
    
    # Reference directory in processing folder
    os.makedirs(os.path.join(config['processing_dir'], 'reference'), exist_ok=True)
    
    # Create merged directory for run tag groups
    os.makedirs(os.path.join(config['processing_dir'], 'merged'), exist_ok=True)
    
    print(f"Created directory structure at {config['processing_dir']} and {config['results_dir']}")

def copy_reference_files(config):
    """Copy reference genome index and GTF files to processing directory"""
    # Get source paths from config
    src_genome_index_base = config["genome_index"]
    src_gtf_file = config["gtf_file"]
    
    # Define target paths
    target_dir = os.path.join(config['processing_dir'], 'reference')
    target_genome_base = os.path.join(target_dir, os.path.basename(src_genome_index_base))
    target_gtf_file = os.path.join(target_dir, os.path.basename(src_gtf_file))
    
    # Store the paths in config for use by other rules
    config["processing_genome_index"] = target_genome_base
    config["processing_gtf_file"] = target_gtf_file
    
    # Copy GTF file if it doesn't exist
    if not os.path.exists(target_gtf_file):
        print(f"Copying GTF file from {src_gtf_file} to {target_gtf_file}")
        try:
            shutil.copy2(src_gtf_file, target_gtf_file)
        except Exception as e:
            print(f"Warning: Failed to copy GTF file: {e}")
    
    # Copy all genome index files
    # BT2 index files have extensions like .1.bt2, .2.bt2, etc.
    index_extensions = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2',
                        '.1.bt2l', '.2.bt2l', '.3.bt2l', '.4.bt2l', '.rev.1.bt2l', '.rev.2.bt2l']
    
    for ext in index_extensions:
        src_file = src_genome_index_base + ext
        target_file = target_genome_base + ext
        
        if os.path.exists(src_file) and not os.path.exists(target_file):
            print(f"Copying genome index file from {src_file} to {target_file}")
            try:
                shutil.copy2(src_file, target_file)
            except Exception as e:
                print(f"Warning: Failed to copy genome index file {src_file}: {e}")

def create_sample_directories(config, sample_name):
    """Create the necessary directory structure for a specific sample"""
    if not sample_name or not isinstance(sample_name, str):
        print(f"Warning: Invalid sample name: {sample_name}, skipping directory creation")
        return
        
    base_dir = os.path.join(config['processing_dir'], sample_name)
    
    # Main directories
    os.makedirs(base_dir, exist_ok=True)
    
    # QC directory
    qc_path = os.path.join(base_dir, 'qc')
    os.makedirs(qc_path, exist_ok=True)
    
    # Logs directory
    logs_dir = os.path.join(base_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    
    # Benchmarks directory
    benchmark_dir = os.path.join(base_dir, "benchmarks")
    os.makedirs(benchmark_dir, exist_ok=True)

def create_run_tag_directories(config, run_tag):
    """Create the necessary directory structure for a run tag"""
    if not run_tag or not isinstance(run_tag, str):
        print(f"Warning: Invalid run tag: {run_tag}, skipping directory creation")
        return
        
    # Create merged directory for BAM merging
    merged_dir = os.path.join(config['processing_dir'], 'merged', run_tag)
    os.makedirs(merged_dir, exist_ok=True)
    os.makedirs(os.path.join(merged_dir, 'logs'), exist_ok=True)
    os.makedirs(os.path.join(merged_dir, 'benchmarks'), exist_ok=True)
    
    # Create main directory for downstream processing
    base_dir = os.path.join(config['processing_dir'], run_tag)
    os.makedirs(base_dir, exist_ok=True)
    
    # QC directory
    qc_path = os.path.join(base_dir, 'qc')
    os.makedirs(qc_path, exist_ok=True)
    os.makedirs(os.path.join(qc_path, 'markduplicates'), exist_ok=True)
    
    # Logs directory
    logs_dir = os.path.join(base_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    
    # Benchmarks directory
    benchmark_dir = os.path.join(base_dir, "benchmarks")
    os.makedirs(benchmark_dir, exist_ok=True)

def load_samples(config):
    """Load sample information from CSV file and organize run tags"""
    try:
        print(f"Loading samples from {config['samples_csv']}")
        
        if not os.path.exists(config['samples_csv']):
            error_msg = f"Sample file {config['samples_csv']} not found!"
            print(f"Error: {error_msg}")
            raise FileNotFoundError(error_msg)
        
        
        samples_df = pd.read_csv(config['samples_csv'])
        
        print(f"Columns in sample file: {', '.join(samples_df.columns)}")
        
        # Validate required columns
        required_columns = ['sample_name', 'read_type', 'source_type', 'file_path_1']
        missing_columns = [col for col in required_columns if col not in samples_df.columns]
        
        if missing_columns:
            error_msg = f"Missing required columns in sample CSV: {', '.join(missing_columns)}"
            print(f"\033[91mError: {error_msg}\033[0m")
            raise ValueError(error_msg)
        
        if samples_df['sample_name'].duplicated().any():
            duplicates = samples_df['sample_name'][samples_df['sample_name'].duplicated()].unique().tolist()
            error_msg = f"Duplicate sample names found in CSV: {', '.join(duplicates)}"
            print(f"\033[91mError: {error_msg}\033[0m")
            raise ValueError(error_msg)        
        
        
        # Validate read_type values
        valid_read_types = ['single', 'paired']
        invalid_read_types = samples_df[~samples_df['read_type'].isin(valid_read_types)]['read_type'].unique()
        
        if len(invalid_read_types) > 0:
            print(f"Error: Invalid read_type values: {', '.join(map(str, invalid_read_types))}. Must be one of: {', '.join(valid_read_types)}")
            raise ValueError(f"Invalid read_type values: {', '.join(map(str, invalid_read_types))}. Must be one of: {', '.join(valid_read_types)}")
        
        # Validate source_type values
        valid_source_types = ['local', 'ftp']
        invalid_source_types = samples_df[~samples_df['source_type'].isin(valid_source_types)]['source_type'].unique()
        
        if len(invalid_source_types) > 0:
            print(f"Error: Invalid source_type values: {', '.join(map(str, invalid_source_types))}. Must be one of: {', '.join(valid_source_types)}")
            raise ValueError(f"Invalid source_type values: {', '.join(map(str, invalid_source_types))}. Must be one of: {', '.join(valid_source_types)}")
        
        # Validate file paths based on read type and source type
        for idx, row in samples_df.iterrows():
            if row['read_type'] == 'paired':
                if pd.isna(row.get('file_path_1', None)) or pd.isna(row.get('file_path_2', None)):
                    print(f"Error: Missing file_path_1 or file_path_2 for paired-end sample {row['sample_name']}")
                    raise ValueError(f"Missing file_path_1 or file_path_2 for paired-end sample {row['sample_name']}")
            else:  # single-end
                if pd.isna(row.get('file_path_1', None)):
                    print(f"Error: Missing file_path_1 for single-end sample {row['sample_name']}")
                    raise ValueError(f"Missing file_path_1 for single-end sample {row['sample_name']}")
        
        # Filter samples if specified
        if config['selected_samples']:
            if isinstance(config['selected_samples'], str):
                # Handle comma-separated string from command line
                selected = config['selected_samples'].split(',')
                samples_df = samples_df[samples_df['sample_name'].isin(selected)]
            else:
                # Handle list from config file
                samples_df = samples_df[samples_df['sample_name'].isin(config['selected_samples'])]
            
            if samples_df.empty:
                print("Error: No samples match the selected_samples filter")
                raise ValueError("No samples match the selected_samples filter")
        
        # Process run tag information if available
        run_tag_samples = {}
        run_tags = []
        effective_samples = []
        
        # Set sample_name as index for easier lookup
        samples_df.set_index('sample_name', inplace=True)
        
        # Check if run_tag column exists
        has_run_tags = 'run_tag' in samples_df.columns
        
        if has_run_tags:
            print("Run tag column found, organizing samples by run tag...")
            
            # Group samples by run tag
            for sample, row in samples_df.iterrows():
                run_tag = row.get('run_tag', '')
                
                # If run tag is empty or NaN, use sample name as run tag
                if pd.isna(run_tag) or run_tag == '':
                    run_tag = sample
                
                # Add run tag to effective samples if not already there
                if run_tag not in effective_samples:
                    effective_samples.append(run_tag)
                    
                # Add run tag to list of run tags if it's different from the sample name
                if run_tag != sample and run_tag not in run_tags:
                    run_tags.append(run_tag)
                
                # Add sample to the run tag group
                if run_tag not in run_tag_samples:
                    run_tag_samples[run_tag] = []
                run_tag_samples[run_tag].append(sample)
            
            # Validate that all samples in a run tag group have the same read type
            for run_tag, samples in run_tag_samples.items():
                if len(samples) > 1:
                    read_types = set(samples_df.loc[samples, 'read_type'])
                    if len(read_types) > 1:
                        print(f"Error: Run tag '{run_tag}' contains samples with different read types: {read_types}")
                        raise ValueError(f"Run tag '{run_tag}' contains samples with different read types: {read_types}")
            
            # Create directories for run tags
            for run_tag in run_tags:
                create_run_tag_directories(config, run_tag)
                print(f"  Run tag: {run_tag} - Contains {len(run_tag_samples[run_tag])} samples: {', '.join(run_tag_samples[run_tag])}")
        else:
            # If no run tag column, each sample is its own effective sample
            effective_samples = list(samples_df.index)
            for sample in effective_samples:
                run_tag_samples[sample] = [sample]
        
        # Create directories for all samples
        print(f"Loaded {len(samples_df)} samples:")
        for sample, row in samples_df.iterrows():
            sample_type = "paired" if row['read_type'] == 'paired' else "single"
            source_type = row['source_type']
            print(f"  {sample}: {sample_type}-end, {source_type} source")
            create_sample_directories(config, sample)
        
        return samples_df, run_tag_samples, run_tags, effective_samples
    
    except Exception as e:
        print(f"Error loading samples: {str(e)}")
        # Create empty DataFrame for graceful failure
        print("Creating empty DataFrame as fallback")
        return pd.DataFrame(columns=['read_type', 'source_type', 'file_path_1', 'file_path_2']).set_index(pd.Index([])), {}, [], []

def validate_genome_index(config):
    """Validate that the reference genome index files exist"""
    index_base = config["genome_index"]
    index_extensions = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
    
    # Check for at least one index file
    found_any = False
    for ext in index_extensions:
        if os.path.exists(index_base + ext):
            found_any = True
            break
    
    if not found_any:
        raise FileNotFoundError(
            f"ERROR: Genome index not found at {index_base}\n"
            f"Please check that the genome_index path in config.yaml is correct.\n"
            f"Expected to find files with extensions: {', '.join(index_extensions)}"
        )
    
    print(f"Found genome index at {index_base}")

def is_sample_completed(config, sample, samples_df):
    """Check if a sample has already been processed"""
    try:
        # Get the effective name for this sample (sample itself or its run tag)
        effective_name = sample
        if 'run_tag' in samples_df.columns:
            run_tag = samples_df.loc[sample, 'run_tag']
            if not pd.isna(run_tag) and run_tag != '':
                effective_name = run_tag
        
        # Check for final results
        final_results_exist = True
        
        # For paired-end samples
        if samples_df.loc[sample, 'read_type'] == 'paired':
            # Check if counts files exist
            paired_counts_all = os.path.join(config['results_dir'], effective_name, f"{effective_name}_counts_paired_all.txt")
            paired_counts_unique = os.path.join(config['results_dir'], effective_name, f"{effective_name}_counts_paired_unique.txt")
            
            # Check if coverage files exist
            all_bw = os.path.join(config['results_dir'], effective_name, f"{effective_name}_all.bw")
            unique_bw = os.path.join(config['results_dir'], effective_name, f"{effective_name}_unique.bw")
            
            # Check if all required output files exist
            final_results_exist = (
                os.path.exists(paired_counts_all) and 
                os.path.exists(paired_counts_unique) and
                os.path.exists(all_bw) and 
                os.path.exists(unique_bw)
            )
        else:  # single-end samples
            # Check if counts files exist
            single_counts_all = os.path.join(config['results_dir'], effective_name, f"{effective_name}_counts_single_all.txt")
            single_counts_unique = os.path.join(config['results_dir'], effective_name, f"{effective_name}_counts_single_unique.txt")
            
            # Check if coverage files exist
            all_bw = os.path.join(config['results_dir'], effective_name, f"{effective_name}_all.bw")
            unique_bw = os.path.join(config['results_dir'], effective_name, f"{effective_name}_unique.bw")
            
            # Check if all required output files exist
            final_results_exist = (
                os.path.exists(single_counts_all) and 
                os.path.exists(single_counts_unique) and
                os.path.exists(all_bw) and 
                os.path.exists(unique_bw)
            )
        
        return final_results_exist
    except Exception as e:
        print(f"Error checking if sample {sample} is completed: {str(e)}")
        return False