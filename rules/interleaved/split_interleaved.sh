#!/bin/bash

# Script to split interleaved FASTQ.gz file into separate R1 and R2 files with modified headers
# Usage: ./split_interleaved_fastq.sh input_interleaved.fastq.gz output_prefix

# Check if correct number of arguments provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_interleaved.fastq.gz> <output_prefix>"
    echo "Example: $0 sample_interleaved.fastq.gz sample"
    echo "Output: sample_1.fq.gz and sample_2.fq.gz"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_PREFIX="$2"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found!"
    exit 1
fi

echo "Processing interleaved FASTQ file: $INPUT_FILE"
echo "Output prefix: $OUTPUT_PREFIX"

# Create temporary directory for intermediate files
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

# Split the interleaved file into separate R1 and R2 files
echo "Splitting interleaved file..."
zcat "$INPUT_FILE" | awk '
BEGIN {
    r1_file = "'$TEMP_DIR'/temp_R1.fastq"
    r2_file = "'$TEMP_DIR'/temp_R2.fastq"
}
{
    # Every 8 lines represents one pair (4 lines each)
    line_in_pair = (NR - 1) % 8
    
    if (line_in_pair < 4) {
        # First read (R1)
        print $0 > r1_file
    } else {
        # Second read (R2)
        print $0 > r2_file
    }
}'

echo "Modifying headers and compressing..."

# Process R1 file: modify headers and compress
awk '{print (NR%4 == 1) ? "@1_" ++i " READ/1": $0}' "$TEMP_DIR/temp_R1.fastq" | gzip > "${OUTPUT_PREFIX}_1.fq.gz"

# Process R2 file: modify headers and compress
awk '{print (NR%4 == 1) ? "@1_" ++i " READ/2": $0}' "$TEMP_DIR/temp_R2.fastq" | gzip > "${OUTPUT_PREFIX}_2.fq.gz"

echo "Done!"
echo "Created files:"
echo "  ${OUTPUT_PREFIX}_1.fq.gz"
echo "  ${OUTPUT_PREFIX}_2.fq.gz"

# Optional: show some statistics
R1_READS=$(zcat "${OUTPUT_PREFIX}_1.fq.gz" | wc -l | awk '{print $1/4}')
R2_READS=$(zcat "${OUTPUT_PREFIX}_2.fq.gz" | wc -l | awk '{print $1/4}')

echo "Statistics:"
echo "  R1 reads: $R1_READS"
echo "  R2 reads: $R2_READS"