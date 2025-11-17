#!/bin/bash

IN_VCF=$1
REF=

# strip file name + name outs
OUT_BASE=$(basename $IN_VCF .vcf.gz)
NORMED_OUT="normed.$OUT_BASE.vcf.gz"
ATOM_OUT="atomed.normed.$OUT_BASE.vcf.gz"
SORTED_OUT="sorted.atomed.normed.$OUT_BASE.vcf.gz"

# normalization
echo "Normalizing vcf..."
bcftools norm -m- -f $REF $IN_VCF -o $NORMED_OUT -O z
echo "Finished normalization"

# atomization
echo "Atomizing vcf..."
bcftools norm -a $NORMED_OUT -o $ATOM_OUT -O z
bcftools index $ATOM_OUT
echo "Finished atomizing"

# sort
echo "Sorting vcf..."
bcftools sort $ATOM_OUT -o $SORTED_OUT -O z
bcftools index $SORTED_OUT
echo "Finished sorting"

# remove intermediates
rm $NORMED_OUT*
rm $ATOM_OUT*