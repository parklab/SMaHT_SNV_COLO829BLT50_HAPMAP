#!/bin/bash

# This script removes duplicate reds from a BAM file in the provided interval, saves into a new BAM file and creates the index

while [ "$1" != "" ]; do
    case $1 in
    --input_bam)
        shift
        input_bam=$1
        ;;
    --processes)
        shift
        processes=$1
        ;;
    --output)
        shift
        output=$1
        ;;
    --reference)
        shift
        reference=$1
        ;;
    --interval_bed)
        shift
        interval_bed=$1
        ;;
    --job)
        shift
        job=$1
        ;;
    esac
    shift
done

echo "JOB $job"
echo "samtools view"

samtools view --region-file $interval_bed -@ $processes -b -F 0x400 $input_bam > ${output}.removed_dups.bam || exit 1

echo "smamtools index"
samtools index ${output}.removed_dups.bam || exit 1