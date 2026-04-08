#!/bin/bash
date
DataPath=$1
BL_bam=$2
T_bam=$3
REF=$4
Interval=$5
Interval_idx=$6
OutPath=$7
IDX=$8
i=$9

echo $DataPath
echo $REF

Manta=/n/data1/hms/dbmi/park/SOFTWARE/Manta/manta-1.6.0.centos6_x86_64/bin
$Manta/configManta.py \
--normalBam $DataPath/$BL_bam \
--tumorBam $T_bam \
--referenceFasta $REF \
--callRegions $Interval \
--runDir $OutPath/$IDX'_'$i 

/n/app/python/2.7.12/bin/python $OutPath/$IDX'_'$i/runWorkflow.py
date
