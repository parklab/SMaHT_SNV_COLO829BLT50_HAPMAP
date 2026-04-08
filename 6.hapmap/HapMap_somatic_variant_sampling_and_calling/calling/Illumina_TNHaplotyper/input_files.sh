#!/bin/bash

# Variables
sample_folder_path=$1
filename_in=$2

# Functions
generate_paths_with_file() {
    local filename="$1"
    input_files=""

    # Loop from 0 to 61
    for i in {0..61}
    do  
        # Print the full path
        input_files+=" ${sample_folder_path}/TNH2_${i}/${filename}"
    done

    echo "$input_files"
}

generate_paths_with_file $filename_in
