#!/bin/bash
log_dir="/n/data1/hms/dbmi/park/jiny/SMaHT/COLO829/1.benchmark/workflow/logs"
mkdir -p "$log_dir"

# Define the base directory containing all subdirectories
base_dir="/n/data1/hms/dbmi/park/jiny/SMaHT/COLO829/Benchmark/Variants/Challenge_2025/AdxBLT/SNV"

# Define the germline VCF file
germline_vcf="/n/data1/hms/dbmi/park/dominika/testing/smaht/haplotyper/negative_control_colo829/colo829bl/output.vcf.gz"

# Define the directory containing the BED files
bed_dir="/n/data1/hms/dbmi/park/jiny/SMaHT/COLO829/1.benchmark/nc_snv"

# Get the list of BED files
bed_files=("$bed_dir"/*.bed)
num_bed_files=${#bed_files[@]}

# Print the number of BED files
echo "Number of BED files: $num_bed_files"

# Loop over all subdirectories in the base directory
for sub_dir in "$base_dir"/*/;
do
    echo "DEBUG: Looking for VCFs in $sub_dir"
    ls "$sub_dir"/*norm.snps.vcf.gz
    for input_vcf in "$sub_dir"*norm.snps.vcf.gz;
    do

        # Define the output file names
        output_germline_vcf="${input_vcf%.vcf}.NC.germline.vcf"
        output_bed_vcf="${input_vcf%.vcf}.NC.nonvariant.vcf"
        echo $input_vcf "submitted"
        # Create a SLURM job script
        job_script=$(mktemp)
        cat <<EOF > "$job_script"
#!/bin/bash
#SBATCH --job-name=intersect_job
#SBATCH -A park_contrib
#SBATCH -p priopark
#SBATCH -w compute-p-17-[174-182]
#SBATCH -N 1
#SBATCH --output=${log_dir}/$(basename ${input_vcf%.vcf.gz}).%j.out
#SBATCH --error=${log_dir}/$(basename ${input_vcf%.vcf.gz}).%.j.err
#SBATCH --time=03:00:00
#SBATCH --mem=2G

# Load the necessary modules
module load bedtools
#bcftools="/n/data1/hms/dbmi/park/jiny/smaht_env/bcftools"

tabix -f ${input_vcf}

# Get the sample name (file name without extension)
sample_name=\$(basename "${input_vcf%.vcf.gz}")

# Create output directory for bcftools isec
isec_output_dir="${sub_dir}\${sample_name}_isec"
mkdir -p "\$isec_output_dir"

# Run bcftools isec
/n/data1/hms/dbmi/park/jiny/smaht_env/bcftools isec -p "\$isec_output_dir" "${input_vcf}" "$germline_vcf"

# Move the generated 0002.vcf to the desired output file
mv "\$isec_output_dir/0002.vcf" "$output_germline_vcf"
echo "$input_vcf intersected with $germline_vcf to make $output_germline_vcf"

# Create a temporary file to store the list of BED files
bed_list_file=\$(mktemp)
ls "$bed_dir"/*_sorted.m.bed > \$bed_list_file

# Initialize a flag for the first BED file
first_bed_file=true
# Loop over all _sorted_merged.bed files and append to the bed output file
while read -r bed_file; 
do
    echo "Processing with BED file: \$bed_file"
    if [ "\$first_bed_file" = true ]; then
        bedtools intersect -a "$input_vcf" -b "\$bed_file" -filenames > "$output_bed_vcf"
        first_bed_file=false
    else
        bedtools intersect -a "$input_vcf" -b "\$bed_file" -filenames >> "$output_bed_vcf"
    fi
done < "\$bed_list_file"
# Print a message indicating that processing of this VCF file is complete
echo "Processing complete for: $input_vcf"

# Clean up the temporary bed list file
rm "\$bed_list_file"
EOF

        # Submit the SLURM job
        sbatch "$job_script"

        # Clean up the job script
        rm "$job_script"
    done
done


