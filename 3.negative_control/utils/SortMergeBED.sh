#!/bin/bash
#SBATCH -c 1
#SBATCH -t 5:00:00
#SBATCH --mem=4G
#SBATCH -p short
#SBATCH -o ./slurm/%j.out
#SBATCH -e ./slurm/%j.out
#SBATCH --array=1-607
module load bedtools

for bed in *_${SLURM_ARRAY_TASK_ID}.bed; do
    bedtools sort -i "$bed" | bedtools merge -i - > merged_sorted/"$bed"
done