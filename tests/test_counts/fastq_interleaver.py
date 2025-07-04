#!/usr/bin/env python3
"""
FASTQ Interleaver Script
Combines two paired FASTQ.gz files into a single interleaved FASTQ file.
"""

import gzip
import argparse
import sys
from pathlib import Path

def read_fastq_record(file_handle):
    """Read a single FASTQ record (4 lines) from file handle."""
    header = file_handle.readline().decode('utf-8').strip()
    if not header:
        return None
    
    sequence = file_handle.readline().decode('utf-8').strip()
    plus = file_handle.readline().decode('utf-8').strip()
    quality = file_handle.readline().decode('utf-8').strip()
    
    return (header, sequence, plus, quality)

def write_fastq_record(file_handle, record):
    """Write a FASTQ record to file handle."""
    header, sequence, plus, quality = record
    file_handle.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")

def interleave_fastq(file1_path, file2_path, output_path):
    """
    Interleave two paired FASTQ.gz files into a single FASTQ file.
    
    Args:
        file1_path: Path to first FASTQ.gz file (R1)
        file2_path: Path to second FASTQ.gz file (R2)
        output_path: Path for output interleaved FASTQ file
    """
    
    # Validate input files exist
    if not Path(file1_path).exists():
        raise FileNotFoundError(f"File not found: {file1_path}")
    if not Path(file2_path).exists():
        raise FileNotFoundError(f"File not found: {file2_path}")
    
    record_count = 0
    
    try:
        with gzip.open(file1_path, 'rb') as f1, \
             gzip.open(file2_path, 'rb') as f2, \
             open(output_path, 'w') as out:
            
            print(f"Processing {file1_path} and {file2_path}...")
            
            while True:
                # Read one record from each file
                record1 = read_fastq_record(f1)
                record2 = read_fastq_record(f2)
                
                # Check if we've reached end of files
                if record1 is None and record2 is None:
                    break
                elif record1 is None or record2 is None:
                    raise ValueError("Input files have different numbers of records")
                
                # Write records in interleaved order (R1, then R2)
                write_fastq_record(out, record1)
                write_fastq_record(out, record2)
                
                record_count += 1
                
                # Progress indicator
                if record_count % 100000 == 0:
                    print(f"Processed {record_count} paired records...")
    
    except Exception as e:
        print(f"Error processing files: {e}", file=sys.stderr)
        # Clean up partial output file
        if Path(output_path).exists():
            Path(output_path).unlink()
        raise
    
    print(f"Successfully interleaved {record_count} paired records")
    print(f"Output written to: {output_path}")

def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Interleave two paired FASTQ.gz files into a single FASTQ file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python fastq_interleaver.py sample_R1.fastq.gz sample_R2.fastq.gz -o sample_interleaved.fastq
  python fastq_interleaver.py reads1.fq.gz reads2.fq.gz -o combined.fastq
        """
    )
    
    parser.add_argument('file1', help='First FASTQ.gz file (R1)')
    parser.add_argument('file2', help='Second FASTQ.gz file (R2)')
    parser.add_argument('-o', '--output', required=True, 
                       help='Output interleaved FASTQ file path')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    try:
        interleave_fastq(args.file1, args.file2, args.output)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()