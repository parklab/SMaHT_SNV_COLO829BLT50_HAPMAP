#!/bin/bash

# Filters every vcf in the source directory by AF intervals of 1% between 0-19, and 20+%

MIN="0"
MAX="0.20"
STEP="0.01"

SOURCE_DIR=./

for file in ${SOURCE_DIR}*vcf.gz; do
    stripped_file="$(basename $file)"
    i=$MIN
    max=0
    mkdir ${stripped_file}_bins_dir
    while [ "$(awk "BEGIN {print ($i < $MAX)}")" -eq 1 ] && [ "$max" != "inf" ]; do
        max=$(awk "BEGIN {print ($i + $STEP)}")

        if awk "BEGIN {exit !($max == $MAX)}"; then
            max="inf"
        fi

        echo "filtering $i - $max" 
        bash filter_af.sh $file af_${i}_${max}_${stripped_file}_vcf $i $max
        bgzip af_${i}_${max}_${stripped_file}_vcf
        bcftools index af_${i}_${max}_${stripped_file}_vcf.gz
        mv af_${i}_${max}_${stripped_file}_vcf.gz* ${stripped_file}_bins_dir

        i=$max
    done
done 
