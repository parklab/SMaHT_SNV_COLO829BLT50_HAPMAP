#!/bin/bash

REF=$1
CONFIDENT_REGIONS=$2

### download  available variant calls ----------------------------------------------------------------------------------------------------------------------

### GIAB cell lines HG002, HG005 
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ChineseTrio/HG005_NA24631_son/latest/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ChineseTrio/HG005_NA24631_son/latest/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi


### pangenome cell lines  HG00438,HG02257,HG02486,HG02622
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc.grch38.vcfbub.a100k.wave.vcf.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc.grch38.vcfbub.a100k.wave.vcf.gz.tbi


### filter calls -------------------------------------------------------------------------------------------------------------------------------------------

# 1.Subset samples in pangenome VCF 
bcftools view -s HG00438,HG02257,HG02486,HG02622 hprc-v1.0-mc.grch38.vcfbub.a100k.wave.vcf.gz -Oz -o hprc-v1.0-mc.grch38.vcfbub.a100k.wave.sample_subset.intermediate.vcf.gz
bcftools index hprc-v1.0-mc.grch38.vcfbub.a100k.wave.sample_subset.intermediate.vcf.gz

# 2. Decompose multi-allelics and MNVs
bcftools norm -a -m -any -Ou HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz| bcftools norm -f $REF -Ou | bcftools sort -Oz -o HG002_GRCh38_1_22_v4.2.1_benchmark.norm.intermediate.vcf.gz
bcftools norm -a -m -any -Ou HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz| bcftools norm -f $REF -Ou | bcftools sort -Oz -o HG005_GRCh38_1_22_v4.2.1_benchmark.norm.intermediate.vcf.gz
bcftools norm -a -m -any -Ou hprc-v1.0-mc.grch38.vcfbub.a100k.wave.sample_subset.intermediate.vcf.gz| bcftools norm -f $REF -Ou | bcftools sort -Oz -o hprc-v1.0-mc.grch38.vcfbub.a100k.wave.sample_subset_norm.intermediate.vcf.gz
bcftools index HG002_GRCh38_1_22_v4.2.1_benchmark.norm.intermediate.vcf.gz
bcftools index HG005_GRCh38_1_22_v4.2.1_benchmark.norm.intermediate.vcf.gz
bcftools index hprc-v1.0-mc.grch38.vcfbub.a100k.wave.sample_subset_norm.intermediate.vcf.gz

# 3. BED region filter 
bcftools view -R $CONFIDENT_REGIONS HG002_GRCh38_1_22_v4.2.1_benchmark.norm.intermediate.vcf.gz -Oz -o HG002_GRCh38_1_22_v4.2.1_benchmark.norm_regionfilter.intermediate.vcf.gz
bcftools view -R $CONFIDENT_REGIONS HG005_GRCh38_1_22_v4.2.1_benchmark.norm.intermediate.vcf.gz -Oz -o HG005_GRCh38_1_22_v4.2.1_benchmark.norm_regionfilter.intermediate.vcf.gz
bcftools view -R $CONFIDENT_REGIONS hprc-v1.0-mc.grch38.vcfbub.a100k.wave.sample_subset_norm.intermediate.vcf.gz -Oz -o hprc-v1.0-mc.grch38.vcfbub.a100k.wave.sample_subset_norm_regionfilter.intermediate.vcf.gz
bcftools index HG002_GRCh38_1_22_v4.2.1_benchmark.norm_regionfilter.intermediate.vcf.gz
bcftools index HG005_GRCh38_1_22_v4.2.1_benchmark.norm_regionfilter.intermediate.vcf.gz
bcftools index hprc-v1.0-mc.grch38.vcfbub.a100k.wave.sample_subset_norm_regionfilter.intermediate.vcf.gz

# 4. Remove indels (keep SNPs only)
bcftools view -v snps HG002_GRCh38_1_22_v4.2.1_benchmark.norm_regionfilter.intermediate.vcf.gz -Oz -o HG002_GRCh38_1_22_v4.2.1_benchmark.norm_regionfilter_SNV.intermediate.vcf.gz
bcftools view -v snps HG005_GRCh38_1_22_v4.2.1_benchmark.norm_regionfilter.intermediate.vcf.gz -Oz -o HG005_GRCh38_1_22_v4.2.1_benchmark.norm_regionfilter_SNV.intermediate.vcf.gz
bcftools view -v snps hprc-v1.0-mc.grch38.vcfbub.a100k.wave.sample_subset_norm_regionfilter.intermediate.vcf.gz -Oz -o hprc-v1.0-mc.grch38.vcfbub.a100k.wave.sample_subset_norm_regionfilter_SNV.intermediate.vcf.gz
bcftools index HG002_GRCh38_1_22_v4.2.1_benchmark.norm_regionfilter_SNV.intermediate.vcf.gz
bcftools index HG005_GRCh38_1_22_v4.2.1_benchmark.norm_regionfilter_SNV.intermediate.vcf.gz
bcftools index hprc-v1.0-mc.grch38.vcfbub.a100k.wave.sample_subset_norm_regionfilter_SNV.intermediate.vcf.gz

# 5. Remove positions with no ALT observations
bcftools view -i 'COUNT(GT~"1")>0' hprc-v1.0-mc.grch38.vcfbub.a100k.wave.sample_subset_norm_regionfilter_SNV.intermediate.vcf.gz -Oz -o hprc-v1.0-mc.grch38.vcfbub.a100k.wave.sample_subset_norm_regionfilter_SNV_wALTallele.intermediate.vcf.gz
bcftools index hprc-v1.0-mc.grch38.vcfbub.a100k.wave.sample_subset_norm_regionfilter_SNV_wALTallele.intermediate.vcf.gz

# 6. Remove any remaining multiallelic sites 
bcftools norm -d all HG002_GRCh38_1_22_v4.2.1_benchmark.norm_regionfilter_SNV.intermediate.vcf.gz -Oz -o HG002_GRCh38_1_22_v4.2.1_benchmark.filtered.vcf.gz
bcftools norm -d all HG005_GRCh38_1_22_v4.2.1_benchmark.norm_regionfilter_SNV.intermediate.vcf.gz -Oz -o HG005_GRCh38_1_22_v4.2.1_benchmark.filtered.vcf.gz
bcftools norm -d all hprc-v1.0-mc.grch38.vcfbub.a100k.wave.sample_subset_norm_regionfilter_SNV_wALTallele.intermediate.vcf.gz -Oz -o hprc-v1.0-mc.grch38.vcfbub.a100k.wave.filtered.vcf.gz
bcftools index HG002_GRCh38_1_22_v4.2.1_benchmark.filtered.vcf.gz
bcftools index HG005_GRCh38_1_22_v4.2.1_benchmark.filtered.vcf.gz
bcftools index hprc-v1.0-mc.grch38.vcfbub.a100k.wave.filtered.vcf.gz

# clean up intermediate files
rm *intermediate.vcf* 



