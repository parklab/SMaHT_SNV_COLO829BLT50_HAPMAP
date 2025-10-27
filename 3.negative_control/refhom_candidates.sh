# bash /n/data1/hms/dbmi/park/jiny/SMaHT/COLO829/2.Github/3.negative_control/refhom_candidates.sh /n/data1/hms/dbmi/park/dominika/testing/smaht/negative_control/WGS/PacBio/COLO829T/Haplotyper_GVCF/multiple_rg/output.vcf /n/data1/hms/dbmi/park/dominika/testing/smaht/negative_control/WGS/PacBio/COLO829BL/multiple_rg/output.vcf

COLO829T_GVCF=$1
COLO829BL_GVCF=$2

module load bedtools


#HaplotypeCaller in GVCF to BED files
python3 /n/data1/hms/dbmi/park/jiny/SMaHT/COLO829/2.Github/3.negative_control/1a.HaplotyperGVCFtoBed_RefHomonly.py $COLO829T_GVCF COLO829T_refhom_candidate_regions.bed
python3 /n/data1/hms/dbmi/park/jiny/SMaHT/COLO829/2.Github/3.negative_control/1a.HaplotyperGVCFtoBed_RefHomonly.py $COLO829BL_GVCF COLO829BL_refhom_candidate_regions.bed

bedtools intersect -a COLO829T_refhom_candidate_regions.bed -b COLO829BL_refhom_candidate_regions.bed > COLO829BL_T_refhom_candidate_regions.intersect.bed

mkdir -p refHomCandidates
#split into 605 jobs
split -n l/605 -d --additional-suffix=.bed COLO829BL_T_refhom_candidate_regions.intersect.bed refHomCandidates/COLO829BL_T_region_homref_cand_

i=0
for f in refHomCandidates/COLO829BL_T_region_homref_cand_*.bed; do
    mv "$f" "refHomCandidates/tmp_${i}.bed"
    i=$((i+1))
done

i=0
for f in refHomCandidates/tmp_*.bed; do
    mv "$f" "refHomCandidates/COLO829BL_T_region_homref_cand_${i}.bed"
    i=$((i+1))
done



