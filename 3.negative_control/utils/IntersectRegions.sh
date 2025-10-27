#!/bin/bash
#SBATCH -c 1
#SBATCH -t 10:00:00
#SBATCH --mem=8G
#SBATCH -p park
#SBATCH -o slurm/intersect_%j.out
#SBATCH -e slurm/intersect_%j.err

module load bedtools

regions=("/n/data1/hms/dbmi/park/jiny/SMaHT/Region_benchmark/Release_Genomic_Stratification_v2/SMaHT_easy_v2.bed.gz" "/n/data1/hms/dbmi/park/jiny/SMaHT/Region_benchmark/Release_Genomic_Stratification_v2/SMaHT_difficult_v2.bed.gz" "/n/data1/hms/dbmi/park/jiny/SMaHT/Region_benchmark/Release_Genomic_Stratification_v2/SMaHT_extreme_v2.bed.gz")
region_names=("/n/data1/hms/dbmi/park/dominika/testing/smaht/COLO829/negative_control/NegativeControlSet/merged_sorted/regions/NegativeControlEasy_v2" "/n/data1/hms/dbmi/park/dominika/testing/smaht/COLO829/negative_control/NegativeControlSet/merged_sorted/regions/NegativeControlDifficult_v2" "/n/data1/hms/dbmi/park/dominika/testing/smaht/COLO829/negative_control/NegativeControlSet/merged_sorted/regions/NegativeControlExtreme_v2")

for r in "${region_names[@]}"; do
    > ${r}.bed
done

for bed in /n/data1/hms/dbmi/park/dominika/testing/smaht/COLO829/negative_control/NegativeControlSet/merged_sorted/Step*.bed; do
    for i in "${!regions[@]}"; do
        bedtools intersect -a "$bed" -b "${regions[$i]}" >> ${region_names[$i]}.bed
    done
done

for r in "${region_names[@]}"; do
    sort -k1,1 -k2,2n ${r}.bed | bedtools merge -i - > ${r}.merged.bed
done
