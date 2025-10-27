#!/bin/bash

REF=$1
OutPath=$2
Output=$3
REGIONS_FILE=$4
I=$5
q=$6
Q=$7
BAM1=$8
BAM2=$9


module load gcc
module load bcftools


time bcftools mpileup \
-f $REF \
$BAM1 $BAM2 \
--threads 1 \
--ignore-RG -d 2000 \
--regions-file $REGIONS_FILE \
-a AD,ADF,ADR,DP \
-Q $Q -q $q > "$OutPath/$Output"_"$I.vcf" && echo "SUCCESS"


bcftools norm -a -m- -f $REF -o "$OutPath/$Output"_"$I.norm.vcf" "$OutPath/$Output"_"$I.vcf"