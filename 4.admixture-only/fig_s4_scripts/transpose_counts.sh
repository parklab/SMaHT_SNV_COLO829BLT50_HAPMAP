#!/bin/bash
# Usage: ./transpose_variants.sh input_file.txt

input_file="$1"

if [ -z "$input_file" ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

# Initialize arrays
declare -A union_pos_1plus union_pos_2plus union_pos_5plus union_pos_8plus union_pos_10plus
declare -A union_homog_neg_1plus union_homog_neg_2plus union_homog_neg_5plus union_homog_neg_8plus union_homog_neg_10plus

# Parse the file
while IFS= read -r line; do
    # Skip empty lines
    [[ -z "$line" ]] && continue

    # Extract degree and control type
    if [[ $line =~ COUNTS\ FOR\ ([0-9]+)\ DEGREE: ]]; then
        degree="${BASH_REMATCH[1]}"
        continue
    fi

    if [[ $line =~ ^[0-9]+\ degree\ ([a-z_]+)\ with\ coverage: ]]; then
        control_type="${BASH_REMATCH[1]}"
        continue
    fi

    # Extract variant counts
    if [[ $line =~ Variants\ validated\ by\ at\ least\ ([0-9]+)\ reads:\ ([0-9]+) ]]; then
        min_reads="${BASH_REMATCH[1]}"
        count="${BASH_REMATCH[2]}"

        # Build array name
                array_name="union_${control_type}_${min_reads}plus"

        # Store the count with degree as index
        case $array_name in
            union_pos_ctrl_1plus) union_pos_1plus[$degree]=$count ;;
            union_pos_ctrl_2plus) union_pos_2plus[$degree]=$count ;;
            union_pos_ctrl_5plus) union_pos_5plus[$degree]=$count ;;
            union_pos_ctrl_8plus) union_pos_8plus[$degree]=$count ;;
            union_pos_ctrl_10plus) union_pos_10plus[$degree]=$count ;;
            union_homog_neg_ctrl_1plus) union_homog_neg_1plus[$degree]=$count ;;
            union_homog_neg_ctrl_2plus) union_homog_neg_2plus[$degree]=$count ;;
            union_homog_neg_ctrl_5plus) union_homog_neg_5plus[$degree]=$count ;;
            union_homog_neg_ctrl_8plus) union_homog_neg_8plus[$degree]=$count ;;
            union_homog_neg_ctrl_10plus) union_homog_neg_10plus[$degree]=$count ;;
        esac
    fi
done < "$input_file"

# Function to print array in the desired format
print_array() {
    local -n arr=$1
    local name=$2
    local values=""

    # Build values string for degrees 1-5
    for degree in 1 2 3 4 5; do
        if [[ -n "${arr[$degree]}" ]]; then
            if [[ -z "$values" ]]; then
                values="${arr[$degree]}"
            else
                values="$values, ${arr[$degree]}"
            fi
        fi
    done

    echo "${name}=[$values]"
}

# Print all arrays
print_array union_pos_1plus "union_pos_1plus"
print_array union_homog_neg_1plus "union_homog_neg_1plus"
print_array union_pos_2plus "union_pos_2plus"
print_array union_homog_neg_2plus "union_homog_neg_2plus"
print_array union_pos_5plus "union_pos_5plus"
print_array union_homog_neg_5plus "union_homog_neg_5plus"
print_array union_pos_8plus "union_pos_8plus"
print_array union_homog_neg_8plus "union_homog_neg_8plus"
print_array union_pos_10plus "union_pos_10plus"
print_array union_homog_neg_10plus "union_homog_neg_10plus"