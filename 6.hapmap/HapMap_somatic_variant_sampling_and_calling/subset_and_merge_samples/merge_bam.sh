#!/bin/bash

IN_DIR=$1
OUT_DIR=$2
OUT_FILE="${OUT_DIR}/merged.bam"

mkdir -p "$OUT_DIR"

### merge bam files in input directory
echo "starting merge at $(date +'%Y-%m-%d %H:%M:%S')"
samtools merge -@ 4 -o "$OUT_FILE" $(find "$IN_DIR" -maxdepth 1 -name "*.bam")
echo "merging bam files completed at $(date +'%Y-%m-%d %H:%M:%S')"

### index merged bam file
samtools index "$OUT_FILE"
echo "index of merged bam completed at $(date +'%Y-%m-%d %H:%M:%S')"

