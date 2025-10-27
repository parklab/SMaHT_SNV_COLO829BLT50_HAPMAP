#!/bin/bash

GCC=$1

# Input files
declare -A dic_BamPath=(["Broad_Illumina"]=/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/HapMap_Mixture/illuminaNovaseq_bulkWgs/seq_data/ ["WashU_Illumina"]=/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/HapMap_Mixture/illuminaNovaseq_bulkWgs/seq_data/ ["BCM_Illumina"]=/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/HapMap_Mixture/illuminaNovaseq_bulkWgs/seq_data/ ["NYGC_Illumina"]=//n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_NYGC/HapMap_Mixture/illuminaNovaseq_bulkWgs/seq_data/ ["BCM_PB"]=/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/HapMap_Mixture/pacbioHifi_bulkWgs/seq_data/ ["Broad_PB"]=/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/HapMap_Mixture/pacbioHifi_bulkWgs/seq_data/ ["WashU_PB"]=/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/HapMap_Mixture/pacbioHifi_bulkWgs/seq_data/ ["UW_PB"]=/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_UW/HapMap_Mixture/pacbioHifi_bulkFiberSeq/seq_data/)


BAM_DIR=${dic_BamPath[$GCC]}

OUTPUT_DIR="/n/data1/hms/dbmi/park/julia/SMaHT/benchmarking/mpileup/pileupVAF/${GCC}"
[ ! -d $OUTPUT_DIR ] && mkdir -p $OUTPUT_DIR

TARGETS="/n/data1/hms/dbmi/park/julia/SMaHT/benchmarking/mosaic_calls/all_TP_and_FP_sorted.vcf.gz"
REF_GENOME="/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh38_smaht/hg38_no_alt.fa"

# Define regions to process
REGION_FILE="/n/data1/hms/dbmi/park/julia/SMaHT/benchmarking/general/ref_genome_intervals/hg38.100intervals_format.bed"
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

