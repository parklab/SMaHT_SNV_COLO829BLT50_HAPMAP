module load gcc
module load samtools

# Input
# Variables
SAMPLE_NAME=$1
OUTPUT_FOLDER=$2
FOLDER=$3
INDEX=$4
# Export
export SENTIEON_LICENSE=license.rc.hms.harvard.edu:8990
export SENTIEON_INSTALL_DIR=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202503.01
export PATH=$PATH:/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202503.01/bin


# Read the input file and store the content in a string
FILES=""
echo $OUTPUT_FOLDER/$SAMPLE_NAME/$FOLDER
cd $OUTPUT_FOLDER/$SAMPLE_NAME/$FOLDER

vcfs=$(bash /n/data1/hms/dbmi/park/dominika/input_files.sh $OUTPUT_FOLDER/$SAMPLE_NAME/$FOLDER output.vcf.gz $INDEX)
echo $vcfs
priors=$(bash /n/data1/hms/dbmi/park/dominika/input_files.sh $OUTPUT_FOLDER/$SAMPLE_NAME/$FOLDER output.priors $INDEX)
contamination=$(bash /n/data1/hms/dbmi/park/dominika/input_files.sh $OUTPUT_FOLDER/$SAMPLE_NAME/$FOLDER output.contamination $INDEX)
# Run reference sample file
 
sentieon_merge_TNfilter.sh $REFERENCE $SAMPLE_NAME $vcfs $priors $contamination