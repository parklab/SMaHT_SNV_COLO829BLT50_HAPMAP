#!/bin/bash

BAM_INPUT=$1
BAM_INPUT_DIR=$(dirname "$BAM_INPUT") 
REF_INPUT_DIR=$2
OUTPUT_DIR=$3
SINGULARITY_CONTAINER_DIR=$4
CONDA_ENV_PATH=$5
THREADS=16
PLATFORM="hifi_revio"

# run the sandbox like this afterward
singularity exec \
  -B ${BAM_INPUT_DIR},${REF_INPUT_DIR},${OUTPUT_DIR} \
  ${SINGULARITY_CONTAINER_DIR}/clairs-to_0.1.0.sif \
  /opt/bin/run_clairs_to \
  --tumor_bam_fn ${BAM_INPUT} \
  --ref_fn ${REF_INPUT_DIR}/hg38_no_alt.fa \
  --threads ${THREADS} \
  --platform ${PLATFORM} \
  --output_dir ${OUTPUT_DIR} \
  --snv_min_af 0 \
  --indel_min_af 0 \
  --min_coverage 1 \
  --qual 4 \
  --conda_prefix ${CONDA_ENV_PATH}
