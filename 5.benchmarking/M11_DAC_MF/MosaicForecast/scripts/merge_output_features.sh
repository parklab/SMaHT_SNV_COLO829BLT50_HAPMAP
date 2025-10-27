#!/bin/bash

# Base directory containing the numbered folders
base_dir=$1
nb=$2
# Output file path
output_file="$base_dir/merged_output.features"

# Initialize the output file
> "$output_file"

# Loop through each folder from 1 to 1155
for i in $(seq 1 $nb); do
    # Path to the current output.features file
    input_file="$base_dir/$i/output.features"

    # Check if the file exists
    if [[ -f "$input_file" ]]; then
        # If the output file is empty, copy the header
        if [[ ! -s "$output_file" ]]; then
            head -n 1 "$input_file" >> "$output_file"
        fi

        # Append the content without the header
        tail -n +2 "$input_file" >> "$output_file"
    else
        # Print the location of the missing file
        echo "Missing file: $input_file"
    fi
done
