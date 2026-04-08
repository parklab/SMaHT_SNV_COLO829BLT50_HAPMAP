#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --account=
#SBATCH --partition=
#SBATCH --ntasks=1
#SBATCH --mem=150G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=1way_lr_pileup

# This script was used to perform pileups in the subject BLT50 bam files.
# This was run five times for each possible degree of intersection.
# E.g. for sites only validated in a single GCC, SITES_VCF pointed to 
# sorted.way_1_combined.vcf.gz - which is the union of variants that 
# were called by RUFUS, Mutect2, or both callers.

POST_PROCESS=/scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/old_ref/analysis/blt50_long_read/pileups/post_process.sh
BAM_PATH="/scratch/ucgd/lustre-core/UCGD_Research/marth_NIH/SMAHT/COLO829/blt50_pacbio_bams/"
bams=(
SMHTCOLO829BLT50-X-X-M45-B001-bcm-SMAFI918XTXQ-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
SMHTCOLO829BLT50-X-X-M45-B001-bcm-SMAFIKDYQOG8-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
SMHTCOLO829BLT50-X-X-M45-B001-bcm-SMAFIQIOUB6M-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
SMHTCOLO829BLT50-X-X-M45-B001-broad-SMAFI5CKEQWJ-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
SMHTCOLO829BLT50-X-X-M45-B001-broad-SMAFI7LW4K57-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
SMHTCOLO829BLT50-X-X-M45-B001-broad-SMAFIGKEKZUJ-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
SMHTCOLO829BLT50-X-X-M45-B001-broad-SMAFIJJFO8TK-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
SMHTCOLO829BLT50-X-X-M45-B001-broad-SMAFIMBXH8YW-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
SMHTCOLO829BLT50-X-X-M45-B001-washu-SMAFITLLY6DL-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
SMHTCOLO829BLT50-X-X-M45-B004-uwsc-SMAFIO3M6CQA-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
)

# These sites vcfs were created with TODO
SITES_VCF=sorted.way_1_combined.vcf.gz
FRAG=$(basename $SITES_VCF .vcf.gz)
SITES_BED="$FRAG.bed"
MAX_DEPTH=100
REFERENCE_FILE=

VCF_NAME=$(basename "$SITES_VCF")

# Extract positions and convert to BED format
bcftools query -f '%CHROM\t%POS\n' "$SITES_VCF" | awk '{print $1"\t"($2-1)"\t"$2}' > "$SITES_BED"

for bam in "${bams[@]}"; do
    NAME="${VCF_NAME}_${bam}"
    bcftools mpileup -Oz -o "${NAME}_max${MAX_DEPTH}.pileup.vcf.gz" --threads 32 -d$MAX_DEPTH -a "FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP" -R $SITES_BED -f $REFERENCE_FILE "${BAM_PATH}${bam}"
    bash $POST_PROCESS "${NAME}_max${MAX_DEPTH}.pileup.vcf.gz"
    mv "${NAME}_max${MAX_DEPTH}.pileup.vcf.gz" raw_pileups
done