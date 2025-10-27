#!/bin/bash

# Input
SAMPLE_NAME=$1 
REF_DIR=$2
REFERENCE=$REF_DIR/hg38_no_alt.fa
OUT_DIR=$3

# Fill in  Sentieon license and directory
#export SENTIEON_LICENSE=
#export SENTIEON_INSTALL_DIR=
export PATH=$PATH:${SENTIEON_INSTALL_DIR}/bin

# The script takes as input a file with a list of BAM files
# Read the input file and store the content in a string
FILES=""

cd ${OUT_DIR}/${SAMPLE_NAME}

vcfs=$(bash input_files.sh ${OUT_DIR}/${SAMPLE_NAME} output.vcf.gz)
echo $vcfs
priors=$(bash input_files.sh ${OUT_DIR}/${SAMPLE_NAME} output.priors)
contamination=$(bash input_files.sh ${OUT_DIR}/${SAMPLE_NAME} output.contamination)

# Run Filter
TNfilter.sh $REFERENCE $SAMPLE_NAME $vcfs $priors $contamination
