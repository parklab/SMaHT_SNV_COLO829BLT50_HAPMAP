#!/bin/bash

# Input
SAMPLE_NAME=$1
REF_DIR=$2
REFERENCE=$REF_DIR/hg38_no_alt.fa
SHARDS=$3
POP_AF=$4
OUT_DIR=$5
input_files=("${@:6}")

# Fill in  Sentieon license and directory
#export SENTIEON_LICENSE=
#export SENTIEON_INSTALL_DIR=
export PATH=$PATH:${SENTIEON_INSTALL_DIR}/bin

INDEX=${SLURM_ARRAY_TASK_ID}
echo $INDEX

# Generate working dir and move
echo ${OUT_DIR}/${SAMPLE_NAME}/TNH2_${INDEX}
mkdir -p ${OUT_DIR}/${SAMPLE_NAME}/TNH2_${INDEX}
cd ${OUT_DIR}/${SAMPLE_NAME}/TNH2_${INDEX}

# Run
TNhaplotyper2.sh $SHARDS $INDEX $REFERENCE $POP_AF $SAMPLE_NAME $input_file
