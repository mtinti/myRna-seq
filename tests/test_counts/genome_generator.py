#!/usr/bin/env python3
import random
import os

def main():
    # Set random seed for reproducibility
    random.seed(42)
    
    # Create a random FASTA file with one entry and 5000 random base pairs
    genome_id = "test_genome"
    sequence_length = 5000
    bases = ['A', 'T', 'G', 'C']
    
    # Generate random sequence
    sequence = ''.join(random.choice(bases) for _ in range(sequence_length))
    
    # Write to FASTA file
    with open("random_genome.fa", "w") as fasta_file:
        fasta_file.write(f">{genome_id}\n")
        
        # Write sequence in chunks of 60 characters for readability
        for i in range(0, len(sequence), 60):
            fasta_file.write(f"{sequence[i:i+60]}\n")
    
    print(f"Created FASTA file: random_genome.fa with {sequence_length} base pairs")
    
    # Define gene regions
    gene_regions = [
        {"name": "gene1", "start": 10, "end": 500},
        {"name": "gene2", "start": 2000, "end": 2500},
        {"name": "gene3", "start": 4000, "end": 4500}
    ]
    
    # Write GTF file
    with open("annotation.gtf", "w") as gtf_file:
        for gene in gene_regions:
            gene_name = gene["name"]
            start = gene["start"]
            end = gene["end"]
            
            # GTF is 1-based, so add 1 to start position
            start_1based = start + 1
            
            # Use the format exactly as requested
            transcript_id = f"{gene_name}.1"
            
            # Transcript entry
            gtf_file.write(f"{genome_id}\tPrediction\ttranscript\t{start_1based}\t{end}\t.\t+\t.\ttranscript_id \"{transcript_id}\"; gene_id \"{gene_name}\";\n")
            
            # Exon entry
            gtf_file.write(f"{genome_id}\tPrediction\texon\t{start_1based}\t{end}\t.\t+\t.\ttranscript_id \"{transcript_id}\"; gene_id \"{gene_name}\";\n")
            
            # CDS entry with gene31 for gene3
            cds_gene_id = "gene31" if gene_name == "gene3" else gene_name
            gtf_file.write(f"{genome_id}\tPrediction\tCDS\t{start_1based}\t{end}\t.\t+\t.\ttranscript_id \"{transcript_id}\"; gene_id \"{cds_gene_id}\";\n")
    
    print(f"Created GTF file: annotation.gtf with {len(gene_regions)} genes")

if __name__ == "__main__":
    main()