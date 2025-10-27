#!/usr/bin/env bash

# *******************************************
# Script to run TNhaplotyper2 on tumor-normal data.
# Generate OrientationBias and ContaminationModel metrics.
# Implemented to run in distributed mode using shards.
# *******************************************

## Functions
generate_random_string() {
    echo "$(cat /dev/urandom | tr -dc 'a-zA-Z' | fold -w 5 | head -n 1)"
}

process_bam_header() {
    # Variable to store replacement @RG args
    replace_rg_args=""

    # Extract @RG lines from the header of the input BAM
    input_bam="$1"
    sample_name="$2"
    rg_lines=$(samtools view -H "$input_bam" | grep "^@RG")

    # Loop through @RG lines and modify ID field
    # Create --replace_rg arguments
    while IFS= read -r rg_line; do
        OLD_SAMPLE_NAME=$(echo $rg_line | grep -o 'SM:[^ ]*' | cut -d':' -f2)
        orig_rg_id=$(echo "$rg_line" | awk -F'\t' '/ID:/ {print $2}' | sed "s/ID://")
        new_rg_id="${orig_rg_id}-$(generate_random_string)"
        new_rg_line=$(echo "$rg_line" | cut -f 2-)
        new_rg_line=$(echo "$new_rg_line" | sed "s/\(SM:\)[^\\t]*/\1${sample_name}/")
        formatted_new_rg_line=$(echo -e "$new_rg_line" | sed 's/\t/\\t/g')
        replace_rg_args+=" --replace_rg ${orig_rg_id}='${formatted_new_rg_line}' "
    done <<< "$rg_lines"

    # Return replacement @RG args
    echo "$replace_rg_args"
}

## Parse command-line arguments
while getopts ":s:i:r:p:t:b:" opt; do
  case $opt in
    s) shards_file="$OPTARG"
    ;;
    i) shard_index="$OPTARG"
    ;;
    r) genome_reference_fasta="$OPTARG"
    ;;
    p) population_allele_frequencies="$OPTARG"
    ;;
    t) sample_name="$OPTARG"
    ;;
    b) sample_bams="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac
done

# Other settings
nt=$SLURM_CPUS_PER_TASK # number of threads to use in computation

# ******************************************
# 1. Process BAM headers and prepare input files
# ******************************************
input_files=""

# Add processed sample BAMs
for sample_bam in ${sample_bams[*]}; do
    replace_args=$(process_bam_header "$sample_bam" "$sample_name")
    input_files+=" $replace_args -i $sample_bam "
done


# ******************************************
# 2. Extract regions from shard
# ******************************************
grep -P "\@${shard_index}\t" "$shards_file" | cut -f 2 > SHARDS_LIST

regions=""
while read -r line; do
    regions+=" --shard $line"
done <SHARDS_LIST


# ******************************************
# 3. Construct TNhaplotyper2 command line
# ******************************************
command="sentieon driver -t $nt -r $genome_reference_fasta $input_files $regions"
command+=" --algo TNhaplotyper2 --tumor_sample $sample_name --germline_vcf $population_allele_frequencies"

echo "PON: " $PON_FLAG
if [ "$PON_FLAG" = "true" ]; then
    command+=" --pon $PON "
fi

command+=" output.vcf.gz"
command+=" --algo OrientationBias --tumor_sample $sample_name output.priors"
command+=" --algo ContaminationModel --tumor_sample $sample_name  -v $population_allele_frequencies output.contamination"

# ******************************************
# 4. Run TNhaplotyper2 command
# ******************************************

eval $command || exit 1
