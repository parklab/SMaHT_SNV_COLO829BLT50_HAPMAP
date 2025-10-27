#!/bin/bash
#$ -cwd

script=/n/data1/hms/dbmi/park/jiny/SMaHT/COLO829/TruthSet/analysis/Indel_pu/annotate_indel_info-wSecondAlt.py
BAM_T=/n/data1/hms/dbmi/park/jiny/SMaHT/COLO829/TruthSet/input/longreads/COLO829T_Hifi.bam
BAM_BL=/n/data1/hms/dbmi/park/jiny/SMaHT/COLO829/TruthSet/input/longreads/COLO829BL_Hifi.bam
DataPath=$1
VCF=$2
python $script \
-i $VCF \
-b $BAM_T

python $script \
-i $VCF \
-b $BAM_BL


#BAM_T=/n/data1/hms/dbmi/park/jiny/SMaHT/COLO829/TruthSet/input/merge_three_bams/COLO829T_Ill_200X.bam
#BAM_BL=/n/data1/hms/dbmi/park/jiny/SMaHT/COLO829/TruthSet/input/merge_three_bams/COLO829BL_Ill_230X.bam
#python $script \
#-i $VCF \
#-b $BAM_T

#python $script \
#-i $VCF \
#-b $BAM_BL
