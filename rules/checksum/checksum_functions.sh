#!/bin/bash
# Common functions for checksum verification

# Function to check MD5 and fail on mismatch
check_md5sum() {
    local file="$1"
    local expected_md5="$2"
    local file_label="$3"
    local flag_file="$4"
    local log_file="$5"
    
    if [ ! -f "$file" ]; then
        echo "$file_label: FAIL - File not found" >> "$flag_file"
        echo "Error: File '$file' not found" >> "$log_file"
        echo "Pipeline will fail for this sample" >> "$log_file"
        exit 1  # Fail the pipeline for this sample
    fi
    
    if [ -z "$expected_md5" ]; then
        echo "$file_label: SKIP - No checksum provided" >> "$flag_file"
        echo "No checksum provided for $file_label, skipping verification" >> "$log_file"
        return 0  # Continue without failing
    fi
    
    local computed_md5=$(md5sum "$file" | awk '{print $1}')
    if [ "$computed_md5" != "$expected_md5" ]; then
        echo "$file_label: FAIL - MD5 mismatch (expected $expected_md5, got $computed_md5)" >> "$flag_file"
        echo "ERROR: MD5 sum mismatch for $file_label" >> "$log_file"
        echo "  Expected: $expected_md5" >> "$log_file"
        echo "  Computed: $computed_md5" >> "$log_file"
        echo "Pipeline will fail for this sample" >> "$log_file"
        exit 1  # Fail the pipeline for this sample
    fi
    
    echo "$file_label: SUCCESS - MD5 matched" >> "$flag_file"
    echo "MD5 sum matched for $file_label" >> "$log_file"
    return 0
}
