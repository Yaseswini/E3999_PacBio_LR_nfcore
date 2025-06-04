#!/usr/bin/env python3
import csv
import os
import sys

# Get input arguments
metadata_path = sys.argv[1]
base_bam_dir = "./data"

# Open the CSV file for reading
with open(metadata_path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    
    # Iterate through each row of the CSV
    for row in reader:
        bam = row['bam_file'].strip()  # BAM file name
        sample = row['sample_id'].strip()  # Sample ID
        condition = row['condition'].strip()  # Condition
        
        # Construct the full BAM file path
        bam_path = os.path.join(base_bam_dir, bam)
        
        # Check if the BAM file exists and print in the correct format
        if os.path.exists(bam_path):
            print(f"{sample}\t{condition}\t{bam_path}")
        else:
            print(f"Warning: BAM file {bam_path} does not exist!")


