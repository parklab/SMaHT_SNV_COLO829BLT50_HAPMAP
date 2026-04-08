#!/bin/bash
# SJG 09Sept2025 for GTM Lab

# Step 0
# This script intersects each of the GCC-respective call sets together - e.g. RUFUS-BCM vs MUTECT-BCM - and then concatenates the output into a union set (1 per GCC - e.g. U-BCM). Each of the five union sets are then intersected with every possible 5-digit bit mask (with a trailing zero for omission from truth set) to create all 31 possible intersection sets. These 31 set intersection counts are output to upset_counts.txt.

# Step 1
# Each intersection is then concatenated by degree into a single variant set per degree. The counts for these concatenated files are output in degree_counts.tsv. These will be used for the experimental and negative control pileups.  

RUFUS_DIR=/scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/old_ref/rufus_runs/
BCM_RUF=${RUFUS_DIR}blt50_1mb_bcm/sorted.atomed.normed.bcm.400x.merged.SMHTCOLO829BLT50.FINAL.no_inherited.SNV_only.hdaf_single_sample.vcf.gz
BROAD_RUF=${RUFUS_DIR}blt50_1mb_broad/sorted.atomed.normed.broad.160x.merged.SMHTCOLO829BLT50.FINAL.no_inherited.SNV_only.hdaf_single_sample.vcf.gz
NYGC_RUF=${RUFUS_DIR}blt50_1mb_nygc/sorted.atomed.normed.SMHTCOLO829BLT50-X-X-M45-A001-nygc-SMAFIC3QAK3U-sentieon_bwamem_202308.01_GRCh38.aligned.sorted.FINAL.no_inherited.SNV_only.hdaf_single_sample.vcf.gz
UW_RUF=${RUFUS_DIR}blt50_1mb_uw/sorted.atomed.normed.uw.350x.merged.SMHTCOLO829BLT50.FINAL.no_inherited.SNV_only.hdaf_single_sample.vcf.gz
WASHU_RUF=${RUFUS_DIR}blt50_1mb_washu/sorted.atomed.normed.washu.500x.merged.SMHTCOLO829BLT50.FINAL.no_inherited.SNV_only.hdaf_single_sample.vcf.gz
RUFUS_SETS=($BCM_RUF $BROAD_RUF $NYGC_RUF $UW_RUF $WASHU_RUF)

MUTECT_DIR=/scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/dac_mutect_calls/pass_only/
BCM_MUT=${MUTECT_DIR}NORMED.PASS_SNPs_ONLY.COLO829BLT-BCM.vcf.gz
BROAD_MUT=${MUTECT_DIR}NORMED.PASS_SNPs_ONLY.COLO829BLT-BROAD.vcf.gz
NYGC_MUT=${MUTECT_DIR}NORMED.PASS_SNPs_ONLY.COLO829BLT-NYGC.vcf.gz
UW_MUT=${MUTECT_DIR}NORMED.PASS_SNPs_ONLY.COLO829BLT-UW.vcf.gz
WASHU_MUT=${MUTECT_DIR}NORMED.PASS_SNPs_ONLY.COLO829BLT-WASHU.vcf.gz
MUTECT_SETS=($BCM_MUT $BROAD_MUT $NYGC_MUT $UW_MUT $WASHU_MUT)

TRUTH_SET=/scratch/ucgd/lustre-core/UCGD_Research/marth_NIH/SMAHT/COLO829/truth_sets/v1_smaht_dac_snv_indel_challenge/truthset_snv_COLO829/SMaHT_COLO829_SNV_truth_set_v1.0.vcf.gz
UPSET_COUNTS="upset_counts.txt"
degree_counts="degree_counts.tsv"

# Function to convert bit mask to corresponding label
bitmask_to_label() {
    local bitmask="$1"
    local labels=("BCM" "Broad" "NYGC" "UW" "WashU" "Truth Set") 
    local result_labels=()
    
    # Check if input is exactly 6 digits
    if [[ ! "$bitmask" =~ ^[01]{6}$ ]]; then
        echo "Error: Input must be exactly 6 binary digits (0 or 1)"
        return 1
    fi
    
    # Process each bit position
    for i in {0..5}; do
        # Extract the bit at position i (0-indexed from left)
        bit="${bitmask:$i:1}"
        
        # If bit is 1, add corresponding label to result
        if [[ "$bit" == "1" ]]; then
            result_labels+=("${labels[$i]}")
        fi
    done
    
    # Join labels with "&" separator
    if [[ ${#result_labels[@]} -eq 0 ]]; then
        echo "No labels (all bits are 0)"
    else
        # Use printf with %s to join array elements with "&"
        local IFS="&"
        echo "${result_labels[*]}"
    fi
}

get_first_one_index() {
    local mask=$1
    local length=${#mask}
    
    for (( i=0; i<length; i++ )); do
        if [[ ${mask:$i:1} == "1" ]]; then
            echo $i
            return
        fi
    done
    
    # If no 1 is found, return -1
    echo -1
}

run_bcftools_isec() {
    local mask=$1
    local combination_type=$2
 
    # Create output directory for this specific mask
    output_dir="isec_output_${combination_type}_${mask}"
    mkdir -p "$output_dir"

    # Run bcftools isec with current bit mask
    bcftools isec -n~"$mask" -c none -p "$output_dir" \
        "bcm/sorted.union.vcf.gz" "broad/sorted.union.vcf.gz" "nygc/sorted.union.vcf.gz" "uw/sorted.union.vcf.gz" "washu/sorted.union.vcf.gz" "$TRUTH_SET"

    # find index of first 1 in mask
    index=$(get_first_one_index "$mask")
    if [[ "$index" == "-1" ]]; then
        echo "ERROR: could not find first non-zero index, exiting..."
        return
    fi
    first_one_vcf="$output_dir/000${index}.vcf"
    bgzip $first_one_vcf
    bcftools index "$first_one_vcf.gz"
  
    # Don't do any AF filtering 
    cp $first_one_vcf.gz $output_dir/combined_vars.vcf.gz
    bcftools index $output_dir/combined_vars.vcf.gz
    
    # Output per combination numbers to a text file for upset plot viz
    label=$(bitmask_to_label $mask $caller)
    echo -en "'$label': " >> $UPSET_COUNTS
    bcftools view -H $output_dir/combined_vars.vcf.gz | wc -l >> $UPSET_COUNTS

    concat_string="$output_dir/combined_vars.vcf.gz"
    echo "$concat_string"
}

# Generates all possible 5-digit masks with trailing 0 for truth-set
generate_bitmasks_array() {
    local k=$1  # combination size
    local n=5  # number of datasets
    local -n result_array=$2  # reference to result array
    
    # Clear the result array
    result_array=()
    
    # Generate all combinations by iterating through all possible bit patterns
    local total_combinations=$((2**n))
    local i=0
    for ((i=1; i<total_combinations; i++)); do  # Start from 1 to skip empty set
        # Count bits in binary representation
        local bit_count=0
        local temp=$i
        while [ $temp -gt 0 ]; do
            if [ $((temp & 1)) -eq 1 ]; then
                ((bit_count++))
            fi
            temp=$((temp >> 1))
        done
        
        # If this combination has exactly k elements
        if [ $bit_count -eq $k ]; then
            # Convert to bitmask string
            local mask=""
            for ((bit=n-1; bit>=0; bit--)); do
                if [ $(((i >> bit) & 1)) -eq 1 ]; then
                    mask="${mask}1"
                else
                    mask="${mask}0"
                fi
            done
            mask="${mask}0"  # Add trailing 0
            result_array+=("$mask")
        fi
    done
}

# Intersect each GCC-respective set
gccs=("bcm" "broad" "nygc" "uw" "washu")
for i in {0..4}; do
    curr_ruf="${RUFUS_SETS[$i]}"
    curr_mut="${MUTECT_SETS[$i]}"
    curr_gcc="${gccs[$i]}"
    bcftools isec -Oz -p $curr_gcc $curr_ruf $curr_mut
    bcftools concat -a -d snps $curr_gcc/0000.vcf.gz $curr_gcc/0001.vcf.gz $curr_gcc/0002.vcf.gz | bcftools sort | bgzip > $curr_gcc/sorted.union.vcf.gz
    bcftools index $curr_gcc/sorted.union.vcf.gz
done
# Do every possible bitmask combination of the unions
for i in {1..4}; do
    echo "Processing $i-way combinations..."
    declare -a masks
    generate_bitmasks_array $i masks
    concat_string=""
    for mask in "${masks[@]}"; do
        concat_curr=$(run_bcftools_isec "$mask" "${i}way")
        concat_string="$concat_string $concat_curr"
    done
    echo "Concatenating all $i-way combinations into single file..."
    bcftools concat -a -d snps -Oz -o "way_${i}_combined.vcf.gz" $concat_string
    echo -en "$i degree\t" >> $degree_counts
    bcftools view -H "way_${i}_combined.vcf.gz" | wc -l >> $degree_counts
    bcftools sort way_${i}_combined.vcf.gz > sorted.way_${i}_combined.vcf
    bgzip sorted.way_${i}_combined.vcf
    bcftools index sorted.way_${i}_combined.vcf.gz
    rm way_${i}_combined.vcf.gz
done

i="5"
echo "Processing $i-way combinations..."
declare -a masks
generate_bitmasks_array $i masks
last_vcf=$(run_bcftools_isec "111110" "${i}way")
cp $last_vcf "way_${i}_combined.vcf.gz"

echo -en "5 degree\t" >> $degree_counts
bcftools view -H "way_${i}_combined.vcf.gz" | wc -l >> $degree_counts
bcftools sort way_${i}_combined.vcf.gz > sorted.way_${i}_combined.vcf
bgzip sorted.way_${i}_combined.vcf
bcftools index sorted.way_${i}_combined.vcf.gz
rm way_${i}_combined.vcf.gz
echo "Finished all possible isec combinations and put counts in $degree_counts"
