"""
Utility functions for sample handling in the RNA-seq pipeline.
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
    
    print(f"Created directory structure at {config['processing_dir']} and {config['results_dir']}")

def copy_reference_files(config):
    """Copy reference FASTA and GTF files to the processing directory."""
    src_fasta = os.path.expanduser(config["reference_fasta"])
    src_gtf_file = os.path.expanduser(config["gtf_file"])

    # Persist expanded paths back into config so downstream rules use the resolved locations
    config["reference_fasta"] = src_fasta
    config["gtf_file"] = src_gtf_file

    target_dir = os.path.join(config["processing_dir"], "reference")

    fasta_basename = os.path.basename(src_fasta)
    fasta_stem, _ = os.path.splitext(fasta_basename)
    if fasta_stem.endswith((".fa", ".fasta", ".fna")):
        fasta_stem = os.path.splitext(fasta_stem)[0]
    target_fasta = os.path.join(target_dir, fasta_basename)
    target_genome_base = os.path.join(target_dir, fasta_stem)
    target_gtf_file = os.path.join(target_dir, os.path.basename(src_gtf_file))

    config["processing_reference_fasta"] = target_fasta
    config["processing_genome_index"] = target_genome_base
    config["processing_gtf_file"] = target_gtf_file

    if not os.path.exists(target_fasta):
        print(f"Copying FASTA file from {src_fasta} to {target_fasta}")
        try:
            shutil.copy2(src_fasta, target_fasta)
        except Exception as e:
            print(f"Warning: Failed to copy FASTA file {src_fasta}: {e}")

    if not os.path.exists(target_gtf_file):
        print(f"Copying GTF file from {src_gtf_file} to {target_gtf_file}")
        try:
            shutil.copy2(src_gtf_file, target_gtf_file)
        except Exception as e:
            print(f"Warning: Failed to copy GTF file: {e}")

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
        valid_read_types = ['single', 'paired', 'nanopore']
        invalid_read_types = samples_df[~samples_df['read_type'].isin(valid_read_types)]['read_type'].unique()
        
        if len(invalid_read_types) > 0:
            print(f"Error: Invalid read_type values: {', '.join(map(str, invalid_read_types))}. Must be one of: {', '.join(valid_read_types)}")
            raise ValueError(f"Invalid read_type values: {', '.join(map(str, invalid_read_types))}. Must be one of: {', '.join(valid_read_types)}")
        
        # Validate source_type values
        valid_source_types = ['local', 'ftp', 'sra_paired']
        invalid_source_types = samples_df[~samples_df['source_type'].isin(valid_source_types)]['source_type'].unique()
        
        if len(invalid_source_types) > 0:
            print(f"Error: Invalid source_type values: {', '.join(map(str, invalid_source_types))}. Must be one of: {', '.join(valid_source_types)}")
            raise ValueError(f"Invalid source_type values: {', '.join(map(str, invalid_source_types))}. Must be one of: {', '.join(valid_source_types)}")
        
        # Validate file paths based on read type and source type
        for idx, row in samples_df.iterrows():
            if row['read_type'] == 'paired':
                if row['source_type'] == 'sra_paired':
                    if pd.isna(row.get('file_path_1', None)):
                        print(f"Error: Missing SRA accession for sample {row['sample_name']}")
                        raise ValueError(f"Missing SRA accession for sample {row['sample_name']}")
                else:
                    if pd.isna(row.get('file_path_1', None)) or pd.isna(row.get('file_path_2', None)):
                        print(f"Error: Missing file_path_1 or file_path_2 for paired-end sample {row['sample_name']}")
                        raise ValueError(f"Missing file_path_1 or file_path_2 for paired-end sample {row['sample_name']}")
            else:  # single-end or nanopore
                if pd.isna(row.get('file_path_1', None)):
                    print(f"Error: Missing file_path_1 for {row['read_type']} sample {row['sample_name']}")
                    raise ValueError(f"Missing file_path_1 for {row['read_type']} sample {row['sample_name']}")
        
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
            
            # Report run tag information
            for run_tag in run_tags:
                print(f"  Run tag: {run_tag} - Contains {len(run_tag_samples[run_tag])} samples: {', '.join(run_tag_samples[run_tag])}")
        else:
            # If no run tag column, each sample is its own effective sample
            effective_samples = list(samples_df.index)
            for sample in effective_samples:
                run_tag_samples[sample] = [sample]
        
        # Print summary of loaded samples
        print(f"Loaded {len(samples_df)} samples:")
        for sample, row in samples_df.iterrows():
            sample_type = row['read_type']
            source_type = row['source_type']
            print(f"  {sample}: {sample_type}, {source_type} source")
        
        return samples_df, run_tag_samples, run_tags, effective_samples
    
    except Exception as e:
        print(f"Error loading samples: {str(e)}")
        # Create empty DataFrame for graceful failure
        print("Creating empty DataFrame as fallback")
        return pd.DataFrame(columns=['read_type', 'source_type', 'file_path_1', 'file_path_2']).set_index(pd.Index([])), {}, [], []

def validate_reference_inputs(config):
    """Validate that the reference FASTA and annotation files exist."""
    fasta_path = os.path.expanduser(config["reference_fasta"])
    config["reference_fasta"] = fasta_path
    if not fasta_path.lower().endswith(".fa"):
        raise ValueError(
            "ERROR: Reference FASTA must have a .fa extension. "
            f"Received: {config['reference_fasta']}"
        )
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(
            "ERROR: Reference FASTA not found at "
            f"{fasta_path}\n"
            "Please set reference_fasta to an existing .fa file (or provide genome_index without the .fa extension)."
        )

    if os.path.getsize(fasta_path) == 0:
        raise FileNotFoundError(
            "ERROR: Reference FASTA is empty at "
            f"{fasta_path}\n"
            "Please provide a FASTA file with sequence data."
        )

    gtf_path = os.path.expanduser(config["gtf_file"])
    config["gtf_file"] = gtf_path
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(
            "ERROR: GTF annotation not found at "
            f"{gtf_path}\n"
            "Please check that the gtf_file path in config.yaml is correct."
        )

    print(f"Found reference FASTA at {fasta_path} and GTF at {gtf_path}")

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
