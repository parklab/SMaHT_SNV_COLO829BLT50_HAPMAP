#!/bin/bash

INPUT_FILE=$1
DOWNSAMPLED_FILE=$2
SUBSAMPLE_SEED=1
FRACTION=$3
START_TASK=${4:-1}  # where to start this process, could be 1, 2 or 3 for downsampling, indexing, stats. defaults to 1


echo "Script started at $(date +'%Y-%m-%d %H:%M:%S')"

if [[ $START_TASK -le 1 ]]; then
  echo "Starting downsample of ${INPUT_FILE} to ${DOWNSAMPLED_FILE} with fraction ${FRACTION} at $(date +'%Y-%m-%d %H:%M:%S')"
  samtools view -bh -@ 4 \
    --subsample ${FRACTION} \
    --subsample-seed ${SUBSAMPLE_SEED} \
    ${INPUT_FILE}  \
    -o ${DOWNSAMPLED_FILE}
  echo "Finished downsample at $(date +'%Y-%m-%d %H:%M:%S')"
fi

if [[ $START_TASK -le 2 ]]; then
  echo "Starting index of ${DOWNSAMPLED_FILE} at $(date +'%Y-%m-%d %H:%M:%S')"
  samtools index ${DOWNSAMPLED_FILE}
  echo "Finished index at $(date +'%Y-%m-%d %H:%M:%S')"
fi

if [[ $START_TASK -le 3 ]]; then
  echo "Starting stats of ${DOWNSAMPLED_FILE} at $(date +'%Y-%m-%d %H:%M:%S')"
  samtools stats ${DOWNSAMPLED_FILE} > ${DOWNSAMPLED_FILE}.stats.txt
  echo "Finished stats at $(date +'%Y-%m-%d %H:%M:%S')"
fi

echo "Script finished at $(date +'%Y-%m-%d %H:%M:%S')"

