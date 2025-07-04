#!/usr/bin/env python3
import random
import os
import re
import argparse
from collections import defaultdict

def parse_fasta(fasta_file):
    """Parse a FASTA file and return a dictionary of sequences."""
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]  # Extract sequence ID
                current_seq = []
            else:
                current_seq.append(line)
    
    if current_id:
        sequences[current_id] = ''.join(current_seq)
    
    return sequences

def parse_gtf(gtf_file):
    """Parse a GTF file and extract CDS regions."""
    cds_regions = defaultdict(list)
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            
            seqid = fields[0]
            start = int(fields[3])  # GTF is 1-based
            end = int(fields[4])
            strand = fields[6]
            
            # Extract gene_id from attributes
            attributes = fields[8]
            gene_id_match = re.search(r'gene_id\s+"([^"]+)"', attributes)
            gene_id = gene_id_match.group(1) if gene_id_match else "unknown"
            
            # Convert to 0-based coordinates for internal processing
            cds_regions[gene_id].append({
                'seqid': seqid,
                'start': start - 1,  # Convert to 0-based
                'end': end,
                'strand': strand,
                'gene_id': gene_id
            })
    
    return cds_regions

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                  'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

def generate_quality_scores(length, min_qual=20, max_qual=40):
    """Generate quality scores in Phred+33 ASCII format."""
    return ''.join(chr(random.randint(min_qual, max_qual) + 33) for _ in range(length))

def simulate_paired_reads(sequence, region, read_length=100, fragment_size=300):
    """Simulate paired-end reads from a genomic region."""
    # Calculate effective region length
    region_length = region['end'] - region['start']
    
    # Ensure the region is large enough for our reads
    if region_length < fragment_size:
        fragment_size = region_length
    
    # Calculate maximum valid starting position for the fragment
    max_start = region['start'] + region_length - fragment_size
    
    # Select random starting position within the region
    fragment_start = random.randint(region['start'], max_start)
    fragment_end = fragment_start + fragment_size
    
    # Extract the fragment from the sequence
    fragment = sequence[fragment_start:fragment_end]
    
    # Generate paired-end reads
    read1 = fragment[:read_length]
    read2_rc = fragment[-read_length:]  # Before reverse complementing
    
    if region['strand'] == '+':
        read2 = reverse_complement(read2_rc)
    else:
        read1, read2 = reverse_complement(read2_rc), read1
    
    # For Picard MarkDuplicates, we need start positions and orientation
    fragment_info = {
        'start': fragment_start,
        'end': fragment_end,
        'strand': region['strand']
    }
    
    return read1, read2, fragment_info

def create_read_id(gene_id, read_num, is_duplicate, duplicate_group=None):
    """Create a read ID with information for tracking duplicates."""
    if is_duplicate and duplicate_group is not None:
        return f"@{gene_id}_read{read_num}_dup{duplicate_group}"
    else:
        return f"@{gene_id}_read{read_num}"

def write_fastq_file(reads, file_path):
    """Write reads to a FASTQ file."""
    with open(file_path, 'w') as f:
        for read_info in reads:
            read_id = read_info['id']
            sequence = read_info['sequence']
            quality = read_info['quality']
            
            f.write(f"{read_id}\n")
            f.write(f"{sequence}\n")
            f.write("+\n")
            f.write(f"{quality}\n")

def main():
    parser = argparse.ArgumentParser(description='Simulate paired-end FASTQ reads from CDS regions.')
    parser.add_argument('--fasta', required=True, help='Path to the input FASTA file')
    parser.add_argument('--gtf', required=True, help='Path to the input GTF file')
    parser.add_argument('--output_prefix', default='simulated', help='Prefix for output FASTQ files')
    parser.add_argument('--read_length', type=int, default=100, help='Length of reads')
    parser.add_argument('--fragment_size', type=int, default=300, help='Size of fragment for paired-end reads')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')
    
    args = parser.parse_args()
    
    # Set random seed for reproducibility
    random.seed(args.seed)
    
    # Parse input files
    sequences = parse_fasta(args.fasta)
    cds_regions = parse_gtf(args.gtf)
    
    # Determine read counts and duplication rates per gene
    read_counts = {
        "gene1": {"total": 10, "duplicates": 5},
        "gene2": {"total": 8, "duplicates": 4},
        "gene3": {"total": 2, "duplicates": 2}
    }
    
    # Prepare output reads
    reads_r1 = []
    reads_r2 = []
    
    # Generate reads for each gene
    for gene_id, regions in cds_regions.items():
        if gene_id not in read_counts:
            continue
        
        # Get read count configuration for this gene
        total_reads = read_counts[gene_id]["total"]
        duplicate_reads = read_counts[gene_id]["duplicates"]
        unique_reads = total_reads - duplicate_reads
        
        # For each region of the gene (typically one CDS per gene in our case)
        for region in regions:
            # Get the sequence for this region
            sequence_id = region['seqid']
            if sequence_id not in sequences:
                print(f"Warning: Sequence {sequence_id} not found in FASTA file")
                continue
            
            sequence = sequences[sequence_id]
            
            # Generate unique reads
            unique_fragments = []
            for i in range(unique_reads):
                read1, read2, fragment_info = simulate_paired_reads(
                    sequence, region, 
                    read_length=args.read_length, 
                    fragment_size=args.fragment_size
                )
                
                read_id = create_read_id(gene_id, i+1, False)
                
                reads_r1.append({
                    'id': f"{read_id}/1",
                    'sequence': read1,
                    'quality': generate_quality_scores(len(read1))
                })
                
                reads_r2.append({
                    'id': f"{read_id}/2",
                    'sequence': read2,
                    'quality': generate_quality_scores(len(read2))
                })
                
                unique_fragments.append(fragment_info)
            
            # Generate duplicate reads (using the same fragment positions as some unique reads)
            if duplicate_reads > 0:
                # If there are no unique reads, generate fragments for duplicates
                if not unique_fragments:
                    for i in range(duplicate_reads):
                        read1, read2, fragment_info = simulate_paired_reads(
                            sequence, region, 
                            read_length=args.read_length, 
                            fragment_size=args.fragment_size
                        )
                        
                        duplicate_group = i + 1
                        read_id = create_read_id(gene_id, unique_reads + i + 1, True, duplicate_group)
                        
                        reads_r1.append({
                            'id': f"{read_id}/1",
                            'sequence': read1,
                            'quality': generate_quality_scores(len(read1))
                        })
                        
                        reads_r2.append({
                            'id': f"{read_id}/2",
                            'sequence': read2,
                            'quality': generate_quality_scores(len(read2))
                        })
                else:
                    # Use existing fragments for duplicates (to ensure they're detected as duplicates)
                    for i in range(duplicate_reads):
                        # Select a fragment to duplicate
                        fragment_idx = i % len(unique_fragments)
                        fragment_info = unique_fragments[fragment_idx]
                        
                        # Extract fragment from sequence using the same coordinates
                        fragment_start = fragment_info['start']
                        fragment_end = fragment_info['end']
                        fragment = sequence[fragment_start:fragment_end]
                        
                        # Create reads with same coordinates but different read IDs
                        read1 = fragment[:args.read_length]
                        read2_rc = fragment[-args.read_length:]
                        
                        if fragment_info['strand'] == '+':
                            read2 = reverse_complement(read2_rc)
                        else:
                            read1, read2 = reverse_complement(read2_rc), read1
                        
                        # Create read ID with duplication info
                        duplicate_group = fragment_idx + 1
                        read_id = create_read_id(gene_id, unique_reads + i + 1, True, duplicate_group)
                        
                        reads_r1.append({
                            'id': f"{read_id}/1",
                            'sequence': read1,
                            'quality': generate_quality_scores(len(read1))
                        })
                        
                        reads_r2.append({
                            'id': f"{read_id}/2",
                            'sequence': read2,
                            'quality': generate_quality_scores(len(read2))
                        })
    
    # Write output FASTQ files
    write_fastq_file(reads_r1, f"{args.output_prefix}_R1.fastq")
    write_fastq_file(reads_r2, f"{args.output_prefix}_R2.fastq")
    
    print(f"Generated {len(reads_r1)} paired-end reads")
    print(f"Output files: {args.output_prefix}_R1.fastq and {args.output_prefix}_R2.fastq")

if __name__ == "__main__":
    main()
