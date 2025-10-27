#!/bin/bash

# Default values
min_dp_inforSNPs=20
n_threads=1


# Function to display usage
usage() {
  echo "Usage: $0 --bam-dir <bam_dir> --output-dir <output_dir> --ref-fasta <ref_fasta> --input-positions <input_positions> [--min-dp <min_dp_inforSNPs>] --umap <umap_mappability> [--threads <n_threads>] " 1>&2
  exit 1
}
# Parse options
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --bam-dir )
        bam_dir=$2
        shift
        shift
        ;;
        --output-dir )
        output_dir=$2
        shift
        shift
        ;;
        --ref-fasta )
        ref_fasta=$2
        shift
        shift
        ;;
        --input-positions_bed )
        input_positions_bed=$2
        shift
        shift
        ;;
        --min-dp )
        min_dp_inforSNPs=$2
        shift
        shift
        ;;
        --umap )
        umap_mappability=$2
        shift
        shift
        ;;
        --sample-name )
        sample_name=$2
        shift
        shift
        ;;
        --threads )
        n_threads=$2
        shift
        shift
        ;;
        --model )
        model=$2
        shift
        shift
        ;;
        * )
        echo "Unknown option: $1"
        usage
        ;;
    esac
done

# Check required parameters
if [ -z "$bam_dir" ] || [ -z "$output_dir" ] || [ -z "$ref_fasta" ] || [ -z "$input_positions_bed" ] || [ -z "$umap_mappability" ]; then
  echo "Missing required arguments" 1>&2
  usage
fi

#input_positions="/tmp/candidates.tsv"

#gunzip -c "$input_positions_bed" | awk -v val="$sample_name" 'BEGIN {OFS="\t"} !/^#/ {gsub(/;/,"\t",$8); print $1, $2-1, $2, $4, $5, val}' > "$input_positions"

#head $input_positions
#echo $samplename
#sleep 120
:<<END
python /usr/local/bin/Phase.py \
$bam_dir \
$output_dir \
$ref_fasta \
$input_positions_bed \
$min_dp_inforSNPs \
$umap_mappability \
$n_threads \
bam
END

#/usr/local/bin/ReadLevel_Features_extraction.py \
#/n/data1/hms/dbmi/park/dominika/testing/smaht/github/MosaicForecast/ReadLevel_Features_extraction.py \
python  /usr/local/bin/ReadLevel_Features_extraction.py \
$input_positions_bed \
$output_dir/output.features \
$bam_dir \
$ref_fasta \
$umap_mappability \
$n_threads \
bam  || exit 1





