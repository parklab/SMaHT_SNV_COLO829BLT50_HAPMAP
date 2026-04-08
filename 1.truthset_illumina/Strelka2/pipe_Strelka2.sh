#!/bin/bash
date
DataPath=$1
BL_bam=$2
T_bam=$3
REF=$4
Interval=$5
OutPath=$6
IDX=$7
i=$8
MantaPath=$9

Strelka2=/n/data1/hms/dbmi/park/SOFTWARE/Strelka/strelka-2.9.2.centos6_x86_64
$Strelka2/bin/configureStrelkaSomaticWorkflow.py \
--normalBam $DataPath/$BL_bam \
--tumorBam $T_bam \
--referenceFasta $REF \
--callRegions $Interval \
--indelCandidates $MantaPath/$IDX'_'$i/results/variants/candidateSmallIndels.vcf.gz \
--runDir $OutPath/$IDX'_'$i

/n/app/python/2.7.12/bin/python $OutPath/$IDX'_'$i/runWorkflow.py -m local -g 10

date

#--indelCandidates $MantaPath/$IDX'_'$i.config/$IDX'_'$i/results/variants/candidateSmallIndels.vcf.gz
