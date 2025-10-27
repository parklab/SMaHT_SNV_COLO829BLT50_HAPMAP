#!/bin/bash

IN_DIR=$1
OUT_DIr=$2

echo "starting merge at $(date +'%Y-%m-%d %H:%M:%S')"

cd $OUT_DIR
samtools merge -@ 4 -o WashU_BCM_NYGC_Illumina_100X_merge.bam $IN_DIR/WashU_Illumina_100X.bam $IN_DIR/BCM_Illumina_100X.bam $IN_DIR/NYGC_Illumina_100X.bam

echo "merge of downsampled 100X Illumina bams completed at $(date +'%Y-%m-%d %H:%M:%S')"
samtools index WashU_BCM_NYGC_Illumina_100X_merge.bam

echo "index of merged bam completed at $(date +'%Y-%m-%d %H:%M:%S')"

