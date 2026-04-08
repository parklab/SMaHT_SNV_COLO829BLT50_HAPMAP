#!/bin/bash

BAM_DIR=$1
OUTPUT_DIR=$2
[ ! -d $OUTPUT_DIR ] && mkdir -p $OUTPUT_DIR

TARGETS=$3 # variant positions
REF_GENOME=$4 #eg "hg38_no_alt.fa"

# Define regions to process
REGION_FILE=$5 # interval bed file of the ref genome
REGION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${REGION_FILE})

# Output file for this region
OUTPUT_FILE="${OUTPUT_DIR}/pileup_${SLURM_ARRAY_TASK_ID}.txt"

echo "bam dir: $BAM_DIR"
echo "bams: ${BAM_DIR}/*.bam"
echo "ref: $REF_GENOME"
echo "regions:$REGION"
echo "target: $TARGETS"
echo "out: $OUTPUT_FILE"

# Run bcftools mpileup
bcftools mpileup \
  -f "${REF_GENOME}" \
  -T "${TARGETS}" \
  -r "${REGION}" \
  -a AD,ADF,ADR,DP \
  -Q 30 \
  -q 30 \
  --ignore-RG -d 3000 \
  -Ov \
  ${BAM_DIR}/*.bam  \
  > "${OUTPUT_FILE}"

