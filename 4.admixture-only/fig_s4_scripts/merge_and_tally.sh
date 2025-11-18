#!/bin/bash
# SJG for GTM Lab 08Sept2025
# This script takes in normalized pileup vcfs from mpileup and merges them into a single merged vcf. 
# Eached merged file is separated by REF entries and ALT entries. The ALT entries only vcf is then 
# intersected with the original RUFUS/MUTECT2 vcf provided to the mpileup sites option (-T). This is 
# because the ALT only vcf will contain ALTs other than the ones that RUFUS called. 
# After we have our intersect-ALT vcf and REF vcf, we then sum the DP and AO fields, and annotate 
# both vcfs with those sums. Finally, we filter the intersect-ALT vcf by N observations (where N is provided as an argument). 
# Thus we are retaining only ALT variants that have at least N counts corroborating the allele.

# Must have pileups in dir named pileup_vcfs and negative control pileups in neg_ctrl_pileup_vcfs

FULL_PATH=/scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/old_ref/analysis/blt50_long_read/dac_mut_ruf_union_pileups/validation_counts/
UNION_PATH=/scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/old_ref/truth_set_v1_overlaps/dac_mut_ruf_unions_no_truth_set/
COUNTS_OUT="$(pwd)/counts.out"
REF=/scratch/ucgd/lustre-labs/marth/resources/references/human/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

# Iterate through every possible degree
for deg in {1..5}; do
    dirs=("pileup_vcfs" "homog_neg_ctrl_pileup_vcfs")
    prefixes=("pos_ctrl" "homog_neg_ctrl")
    rufus_prefixes=("" "")
    for i in 0 1; do
        dir="${dirs[i]}"
        prefix="${prefixes[i]}"
        rufus_prefix="${rufus_prefixes[i]}"

        echo "COUNTS FOR $deg DEGREE:" >> $COUNTS_OUT

        cd ${FULL_PATH}way_$deg/$dir

        # Add each pileup vcf to list
        for vcf in *vcf.gz; do
            echo "$vcf" >> pileup_list.txt
        done

        # Merge pileup vcfs into one, keeping all alt alleles in a single entry
        echo "Merging normed vcfs..."
        bcftools merge -l pileup_list.txt -m none --force-samples -Oz -o merged_${deg}_way.atomed.normed.vcf.gz
        bcftools index merged_${deg}_way.atomed.normed.vcf.gz

        # Left-align and combine multi-allelics to get site counts with coverage
        bcftools norm -f $REF -m -any merged_${deg}_way.atomed.normed.vcf.gz | \
        bcftools norm -m +any -Oz -o multi_allel_comb_merged_${deg}_way.atomed.normed.vcf.gz
        echo -en "$deg degree $prefix with coverage: " >> $COUNTS_OUT
        bcftools view -H multi_allel_comb_merged_${deg}_way.atomed.normed.vcf.gz | wc -l >> $COUNTS_OUT

        # Already atomized all files post pileup 
        atomed_merged_out="merged_${deg}_way.atomed.normed.vcf.gz"

        # Separate ALT and REF entries
        echo "Separating ALTs and REFs..."
        bcftools view -e 'ALT="<*>"' $atomed_merged_out -Oz -o "${prefix}.ALT_only.merged_${deg}_way.atomed.normed.vcf.gz"
        bcftools view -i 'ALT="<*>"' $atomed_merged_out -Oz -o "${prefix}.REF_only.merged_${deg}_way.atomed.normed.vcf.gz"
        bcftools index "${prefix}.ALT_only.merged_${deg}_way.atomed.normed.vcf.gz"
        bcftools index "${prefix}.REF_only.merged_${deg}_way.atomed.normed.vcf.gz"

        # Intersect ALT vcf with RUFUS/Mutect2 unioned vcf
        echo "Intersecting pileup ALT vcf with unioned RUFUS/Mutect2 vcf..."
        mkdir -p alt_isecs
        UNION_VCF="${UNION_PATH}${rufus_prefix}sorted.way_${deg}_combined.vcf.gz"
        bcftools isec -c none -p alt_isecs "${prefix}.ALT_only.merged_${deg}_way.atomed.normed.vcf.gz" $UNION_VCF
        cp alt_isecs/0002.vcf "${prefix}.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf"
        bgzip "${prefix}.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf"
        bcftools index "${prefix}.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf.gz"

        # Move product to above directory
        mv ./${prefix}.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf.gz* ..
        mv ./${prefix}.REF_only.merged_${deg}_way.atomed.normed.vcf.gz* ..

        # Sum ALT AD field across samples and annotate
        echo "Summing ALT observations..." 
        cd ..
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${prefix}.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf.gz | \
                awk -F'\t' '{
            sum = 0
            for(i=5; i<=NF; i++) {
                    if($i != "." && $i != "") {
                        split($i, ad, ",")
                    if(length(ad) >= 2 && ad[2] != "") {
                        sum += ad[2]
                    }
                }
            }
            print $1 "\t" $2 "\t" $3 "\t" $4 "\t" sum
        }' > alt_ad_annotation.txt

        bgzip -f alt_ad_annotation.txt
        tabix -s1 -b2 -e2 alt_ad_annotation.txt.gz

        # Add the header and annotate
        echo '##INFO=<ID=ALT_AD_SUM,Number=1,Type=Integer,Description="Sum of alternate allele depths across samples">' > header.txt

        bcftools annotate \
            -a alt_ad_annotation.txt.gz \
            -h header.txt \
            -c CHROM,POS,REF,ALT,INFO/ALT_AD_SUM \
            ${prefix}.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf.gz \
            -Oz -o ${prefix}.summed.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf.gz

        echo "Filtering by >=1 ALT observations..." 
        bcftools view -i "ALT_AD_SUM >= 1" -Oz -o "${prefix}.1_AD_filtered.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf.gz" ${prefix}.summed.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf.gz

        echo "Filtering by >=2 ALT observations..." 
        bcftools view -i "ALT_AD_SUM >= 2" -Oz -o "${prefix}.2_AD_filtered.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf.gz" ${prefix}.summed.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf.gz

        echo "Filtering by >=5 ALT observations..." 
        bcftools view -i "ALT_AD_SUM >= 5" -Oz -o "${prefix}.5_AD_filtered.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf.gz" ${prefix}.summed.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf.gz

        echo "Filtering by >=8 ALT observations..." 
        bcftools view -i "ALT_AD_SUM >= 8" -Oz -o "${prefix}.8_AD_filtered.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf.gz" ${prefix}.summed.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf.gz

        bcftools view -i "ALT_AD_SUM >= 10" -Oz -o "${prefix}.10_AD_filtered.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf.gz" ${prefix}.summed.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf.gz
        # Tallying counts
        counts=(1 2 5 8 10)
        for count in "${counts[@]}"; do
            echo -en "Variants validated by at least $count reads: " >> $COUNTS_OUT
            bcftools view -H "${prefix}.${count}_AD_filtered.UNION_ALT_only.merged_${deg}_way.atomed.normed.vcf.gz" | wc -l >> $COUNTS_OUT
        done
        echo "" >> $COUNTS_OUT
    done
done
     