#!/bin/bash

MIN=0
MAX="0.20"
STEP="0.01"

SOURCE_DIR=./af_added.masked_blt50.2read_cutoff.way_2.vcf.gz_bins_dir/
OUT_FILE="admixture_vars_by_degree_af_counts.txt"

echo -e "file\tmin\tmax\tbin_count\tlabel" > $OUT_FILE

i=$MIN
max=0
# bash can't do fp arithmetic in for loops...
while [ "$(awk "BEGIN {print ($i < $MAX)}")" -eq 1 ] && [ "$max" != "inf" ]; do
    min=$i
    max=$(awk "BEGIN {print ($i + $STEP)}")

    if awk "BEGIN {exit !($max == $MAX)}"; then
        max="inf"
    fi
    
    blt50_vcf=$(echo "${SOURCE_DIR}af_${min}_${max}*vcf.gz")
    blt50_count=$(bcftools view -H $blt50_vcf | wc -l)
    min_perc=$(awk "BEGIN {print ($min * 100)}")
    max_perc=$(awk "BEGIN {print ($max * 100)}")
    label="${min_perc}-${max_perc}%"
    file=$(basename $blt50_vcf)
    echo -e "$file\t$min\t$max\t$blt50_count\t$label" >> $OUT_FILE
   
    # increment loop
    i=$max
done
