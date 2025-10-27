#!/bin/bash
#SBATCH -c 16
#SBATCH -t 10:00:00                     # D-HH:MM
#SBATCH -p priopark
#SBATCH -A park_contrib
#SBATCH --mem=32G
#SBATCH -o /n/scratch/users/d/dm334/slurm/%a.log
#SBATCH -e /n/scratch/users/d/dm334/slurm/%a.log
#SBATCH --array=1-605

# This script runs Sentieon Haplotyper (GATK HaplotypeCaller) in GVCF mode for a single sample
# This script was run for Illumina COLO829BL and COLO829T to get homozygous-reference sites for negative control and variant sites for further analysis


while [ "$1" != "" ]; do
    case $1 in
    --input_bam)
        shift
        input_bam=$1
        ;;
    --processes)
        shift
        processes=$1
        ;;
    --output)
        shift
        output=$1
        ;;
    --reference)
        shift
        reference=$1
        ;;
    --known_snps)
        shift
        known_snps=$1
        ;;
    --shard)
        shift
        shards_file=$1
        ;;
    --sample_name)
        shift
        sample_name=$1
        ;;
    esac
    shift
done


job=${SLURM_ARRAY_TASK_ID}

# ******************************************
# 2. Create shards
# ******************************************
SHARDS_LIST="SHARDS_LIST_"${job}
grep -P "\@${job}\t" $shards_file | cut -f 2 > $SHARDS_LIST

shards=""

# Reading shards
while read -r line;
  do
    shards+=" --shard $line"
  done <$SHARDS_LIST


echo "JOB $job"

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

input_files=""

replace_args=$(process_bam_header $input_bam $sample_name)
input_files+=" $replace_args -i $input_bam "



echo "Sentieon Haplotyper"
command="/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202503.01/bin/sentieon driver -t 16 $shards -r $reference $input_files --algo Haplotyper -d $known_snps --emit_mode gvcf ${output}_${job}.vcf.gz"


eval $command || exit 1


echo "SUCCESS $job"