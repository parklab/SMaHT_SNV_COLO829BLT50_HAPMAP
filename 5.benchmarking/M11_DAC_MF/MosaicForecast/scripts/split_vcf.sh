#!/bin/bash

# Define the input VCF file
input_vcf=$1
output_dir=$2
mkdir -p $output_dir
# Extract the header lines
grep '^#' "$input_vcf" > $output_dir/header.vcf

# Split the VCF file into chunks of 100 lines each, excluding the header
grep -v '^#' "$input_vcf" | split -l 1000 - $output_dir/temp_part_

# Rename the split files to have the .vcf extension and prepend the header
counter=1
for file in $output_dir/temp_part_*; do
    new_file="${output_dir}/part_${counter}.vcf"
    mv "$file" "$new_file"
    cat $output_dir/header.vcf "$new_file" > $output_dir/temp.vcf && mv $output_dir/temp.vcf "$new_file"
    counter=$((counter + 1))
    bgzip -f $new_file
    tabix ${new_file}.gz
done

# Clean up the temporary header file
rm $output_dir/header.vcf
echo $counter > $output_dir/JOBS
echo "VCF file has been split into multiple files, each containing 100 lines plus header."