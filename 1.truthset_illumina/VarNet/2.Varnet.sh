#!/bin/bash
#SBATCH -A park_contrib
#SBATCH -p park 
#SBATCH -t 88:00:00
#SBATCH --mem 64G
#SBATCH -c 16
#SBATCH -e /n/data1/hms/dbmi/park/dominika/testing/smaht/varnet/logs/Varnet_%j.err
#SBATCH -o /n/data1/hms/dbmi/park/dominika/testing/smaht/varnet/logs/Varnet_%j.log
#SBATCH --mail-type=END


#This scripts runs VarNet filter.py and predict.py

show_help() {
    echo "Usage: $0 -n <normal_bam> -t <tumor_bam> -i <sample_id> -s <sample_name> -o <output_dir> -p <interval_path>"
    echo "    -n  Path to the normal BAM file"
    echo "    -t  Path to the tumor BAM file"
    echo "    -s  Sample name"
    echo "    -o  Output directory"
    echo "    -p  Path to the intervals directory"
    exit 1
}

# Parsing command-line arguments
while getopts "n:t:i:s:o:p:h" opt; do
    case $opt in
        n) normal_bam="$OPTARG" ;;
        t) tumor_bam="$OPTARG" ;;
        s) sample_name="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        p) interval_path="$OPTARG" ;;
        h) show_help ;;
        *) show_help ;;
    esac
done

ERROR_LOG="/n/data1/hms/dbmi/park/dominika/testing/smaht/varnet/logs/Varnet_${SLURM_JOB_ID}.err"
STDOUT_LOG="/n/data1/hms/dbmi/park/dominika/testing/smaht/varnet/logs/Varnet_${SLURM_JOB_ID}.log"

cd $ouput_dir 

ln -s -f $ERROR_LOG errorLogsBQSR.err
ln -s -f $STDOUT_LOG stdoutLogsBQSR.out

# Check if all required arguments are provided
if [ -z "$normal_bam" ] || [ -z "$tumor_bam" ] || [ -z "$sample_name" ] || [ -z "$output_dir" ] || [ -z "$interval_path" ]; then
    echo "Error: Missing required arguments"
    show_help
fi

normal_path=$(dirname $normal_bam)
echo $normal_path
tumor_path=$(dirname $tumor_bam)
echo $tumor_path
mkdir -p $output_dir

# Run Varnet filter step
singularity exec \
--bind /n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt:/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt \
--bind /n/data1/hms/dbmi/park/dominika/intervals:/n/data1/hms/dbmi/park/dominika/intervals \
--bind /n/data1/hms/dbmi/park-smaht_dac/DATA/:/n/data1/hms/dbmi/park-smaht_dac/DATA/ \
--bind ${output_dir}:${output_dir} \
--bind ${normal_path}:${normal_path} \
--bind ${tumor_path}:${tumor_path} \
/n/app/singularity/containers/park-smaht_dac/varnet_1.1.0.sif \
python3 /varnet/filter.py \
--sample_name $sample_name \
--normal_bam $normal_bam  \
--tumor_bam $tumor_bam \
--processes 16 \
--output_dir ${output_dir} \
--reference /n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa \
--region_bed ${interval_path}

# Run Varnet predict step
singularity exec \
--bind /n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt:/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt \
--bind ${output_dir}:${output_dir} \
--bind ${output_dir}:${output_dir} \
--bind ${normal_path}:${normal_path} \
--bind ${tumor_path}:${tumor_path} \
/n/app/singularity/containers/park-smaht_dac/varnet_1.1.0.sif \
python3 /varnet/predict.py \
--sample_name $sample_name \
--normal_bam $normal_bam \
--tumor_bam $tumor_bam \
--processes 16 \
--output_dir ${output_dir}/ \
--reference /n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
