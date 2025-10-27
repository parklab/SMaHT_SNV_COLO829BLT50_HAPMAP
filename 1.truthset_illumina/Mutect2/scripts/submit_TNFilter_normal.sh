#!/bin/bash
#------


module load gcc
module load samtools

# Input
# Variables
FOLDER=$1
SAMPLE_NAME=$2
NORMAL_NAME=$3
FOLDER_OUTPUT=$4
SHARDS_LEN=$5

# Export
export SENTIEON_LICENSE=license.rc.hms.harvard.edu:8990
export SENTIEON_INSTALL_DIR=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202503.01
export PATH=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202503.01/bin:$PATH


# The script take as input a file with a list of BAM files
# Read the input file and store the content in a string
FILES=""
cd $FOLDER_OUTPUT/$FOLDER

vcfs=$(bash /n/data1/hms/dbmi/park/dominika/input_files.sh $FOLDER_OUTPUT/$FOLDER output.vcf.gz $SHARDS_LEN)
priors=$(bash /n/data1/hms/dbmi/park/dominika/input_files.sh $FOLDER_OUTPUT/$FOLDER output.priors $SHARDS_LEN)
contamination=$(bash /n/data1/hms/dbmi/park/dominika/input_files.sh $FOLDER_OUTPUT/$FOLDER output.contamination $SHARDS_LEN)

# Run reference sample files 
bash /n/data1/hms/dbmi/park/dominika/testing/smaht/tissues/scripts/sentieon_merge_TNfilter_normal.sh $REFERENCE $SAMPLE_NAME $NORMAL_NAME $vcfs $priors $contamination