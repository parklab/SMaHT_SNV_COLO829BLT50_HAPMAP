#!/bin/bash

### merge and norm pileup results

PILEUP_DIR=$1
cd $PILEUP_DIR

PREFIX=$2
SUFFIX=".out"
OUTPUT_FILE="/n/data1/hms/dbmi/park/julia/SMaHT/benchmarking/mpileup/pileupVAF/pileup_all_calls_for_pileupVCFcalc.vcf"
OUTPUT_FILE_NORM="/n/data1/hms/dbmi/park/julia/SMaHT/benchmarking/mpileup/pileupVAF/pileup_all_calls_for_pileupVCFcalc.norm.vcf"

# Find all files matching the pattern and sort them numerically
FILES=$(ls ${PREFIX}*${SUFFIX} 2>/dev/null | sort -V)

# Check if any files were found
if [[ -z "$FILES" ]]; then
    echo "No files found matching the pattern ${PREFIX}*${SUFFIX}."
    exit 1
fi

# Merge the VCF files using bcftools
bcftools concat -o $OUTPUT_FILE -O v $FILES

# Check if the merge was successful
if [[ $? -eq 0 ]]; then
    echo "Files merged successfully into $OUTPUT_FILE."
else
    echo "Error: Failed to merge files."
    exit 1
fi

# Norm the variants using bcftools
bcftools norm -m-any $OUTPUT_FILE > $OUTPUT_FILE_NORM
