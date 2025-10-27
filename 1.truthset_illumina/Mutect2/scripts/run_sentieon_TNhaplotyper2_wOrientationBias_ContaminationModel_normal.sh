#!/bin/bash
#------

module load gcc
module load samtools

# Export
export SENTIEON_LICENSE=license.rc.hms.harvard.edu:8990
export SENTIEON_INSTALL_DIR=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202503.01
export PATH=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202503.01/bin:$PATH

# Command-line arguments
usage() {
    echo "Usage: $0 -s SAMPLE_NAME -t 'SAMPLE_BAM1 SAMPLE_BAM2 ...' -n NORMAL_NAME -m 'NORMAL_BAM1 NORMAL_BAM2 ...' -f FOLDER -i INDEDX"
    exit 1
}

while getopts ":s:t:n:m:f:i:o:" opt; do
    case $opt in
        s) SAMPLE_NAME="$OPTARG" ;;
        t) SAMPLE_BAMS=($OPTARG) ;;
        n) NORMAL_NAME="$OPTARG" ;;
        m) NORMAL_BAMS=($OPTARG) ;;
        f) FOLDER="$OPTARG" ;;
        i) INDEX="$OPTARG" ;;
        o) OUTPUT_FOLDER="$OPTARG" ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [ -z "$SAMPLE_NAME" ] || [ -z "$SAMPLE_BAMS" ] || [ -z "$NORMAL_NAME" ] || [ -z "$NORMAL_BAMS" ] || [ -z "$FOLDER" ] || [ -z "$INDEX" ]; then
    usage
fi

echo $FOLDER


echo "Shard Index: $INDEX"

# PATHS
WORKING_DIR=${OUTPUT_FOLDER}/${SAMPLE_NAME}/${FOLDER}/TNH2_${INDEX}

# Create working directory and move there
echo "Creating and moving to directory: $WORKING_DIR"
mkdir -p "$WORKING_DIR"
cd "$WORKING_DIR" || exit 1

# Run the TNhaplotyper2 script
echo "Running TNhaplotyper2..."
/n/data1/hms/dbmi/park/dominika/testing/smaht/tissues/scripts/sentieon_TNhaplotyper2_wOrientationBias_ContaminationModel_normal.sh \
-s "$SHARDS" -i "$INDEX" -r "$REFERENCE" -p "$POP_AF" -t "$SAMPLE_NAME" -b "${SAMPLE_BAMS[*]}" -n "$NORMAL_NAME" -m "${NORMAL_BAMS[*]}"
