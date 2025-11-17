#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=
#SBATCH --partition=
#SBATCH --ntasks=1
#SBATCH --mem=150G
#SBATCH --cpus-per-task=16 
#SBATCH --job-name=1way_neg_lr_pileup

# This script was used to perform pileups in unrelated benchmarking files that
# utilized the same sequencing and processing pipeline as COLOBLT50 samples,
# to serve as a negative control.
# This was run five times for each possible degree of intersection,
# just like in pileup.slurm.

POST_PROCESS=post_process.sh
BAM_PATH="/scratch/ucgd/lustre-core/UCGD_Research/marth_NIH/SMAHT/benchmarking_tissue_homogenates/DonorST001/liver_1A/pacbio/"
bams=(
ST001-1A-XX-M22-B001-bcm-SMAFI1DZ62NM-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
ST001-1A-XX-M22-B001-bcm-SMAFI9VER4C3-pbmm2_1.13.0_GRCh38.aligned.sorted.bam 
ST001-1A-XX-M22-B001-bcm-SMAFIYT6BC6F-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
ST001-1A-XX-M22-B001-broad-SMAFIFXDVRTL-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
ST001-1A-XX-M22-B001-washu-SMAFIQB3OD1L-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
)

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