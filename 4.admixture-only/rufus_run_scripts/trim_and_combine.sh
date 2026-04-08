#!/bin/bash

#SBATCH --time=0-01:00:00                                                                                                   
#SBATCH --account=marth-rw                                                                                                        
#SBATCH --partition=marth-shared-rw                                                                                              
#SBATCH -o /scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/slurm_out/%j/out.out
#SBATCH -e /scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/slurm_out/%j/error.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=u0746015@utah.edu

# This script concatenates individual vcfs produced by piecemeal RUFUS runs,
# into a single combined vcf. Only variants from within the region targeted by the run
# are kept (at this time, RUFUS may make erroneous calls out of the given range). 
# The final vcf contains the header from the first chr1_1_{CHUNK_SIZE} final 
# individual vcf. Note, the chr and chr length arrays, as well as the chunk size
# must be equal to those run for the piecemeal run.
# Compresses and indexes the final file.

module load bcftools
module load htslib

# Files and dirs for runscript
SOURCE_DIR=/scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/rufus_runs/blt50_1mb_bcm/
RUFUS_VCF=bcm.400x.merged.SMHTCOLO829BLT50.sorted.bam.generator.V2.overlap.hashcount.fastq.bam.FINAL.vcf.gz
TRIMMED_VCF=trimmed.final.vcf

RUFUS_PRE_VCF=Intermediates/bcm.400x.merged.SMHTCOLO829BLT50.sorted.bam.generator.V2.overlap.hashcount.fastq.bam.sorted.vcf
TRIMMED_PRE_VCF=Intermediates/trimmed.prefiltered.vcf

COMBINED_VCF="combined.final.vcf"
COMBINED_PRE_VCF="combined.prefiltered.vcf"

COMBINED_HEADER="/scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/headers/combined.header"
COMBINED_PRE_HEADER="/scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/headers/combined.pre.header"

NUM_CHRS=24
CHUNK_SIZE=1000000

CHRS=(
"1"  
"2" 
"3" 
"4" 
"5" 
"6" 
"7" 
"8" 
"9" 
"10" 
"11" 
"12" 
"13" 
"14" 
"15"
"16" 
"17" 
"18" 
"19" 
"20" 
"21" 
"22" 
"X" 
"Y" 
)

# GRCh38
CHR_LENGTHS=(
248956422
242193529
198295559
190214555
181538259
170805979
159345973
145138636
138394717
133797422
135086622
133275309
114364328
107043718
101991189
90338345
83257441
80373285
58617616
64444167
46709983
50818468
156040895
57227415
)

cd $SOURCE_DIR

# Initialize combined vcf files
cat "$COMBINED_HEADER" > $COMBINED_VCF
cat "$COMBINED_PRE_HEADER" > $COMBINED_PRE_VCF

n=${#CHRS[@]}
for (( i = 0; i < n; i++ ))
do
    curr_len="${CHR_LENGTHS[i]}"
    curr_chr="${CHRS[i]}"
    ((chunks=curr_len/CHUNK_SIZE+1))

    start_coord=1
    for (( c = 1; c <= chunks; c++ ))
    do
        ((end_coord=start_coord+CHUNK_SIZE-1))
        if [ $end_coord -gt $curr_len ]; then
            echo "passed end coordinate, setting to $curr_len"
            end_coord=$curr_len
        fi

        # Get region dir

        if [[ -f "$TARGET_DIR${RUFUS_VCF}" ]]; then

            # Write out trimmed region to final vcf
            bcftools view -H -r "chr${curr_chr}:${start_coord}-${end_coord}" "${TARGET_DIR}${RUFUS_VCF}" >> $COMBINED_VCF

            # Write out trimmed region to prefiltered vcf
            if [[ -f "${TARGET_DIR}${RUFUS_PRE_VCF}" ]]; then
                bgzip "${TARGET_DIR}${RUFUS_PRE_VCF}"
                bcftools index -t "${TARGET_DIR}${RUFUS_PRE_VCF}.gz"
            fi
            bcftools view -H -r "chr${curr_chr}:${start_coord}-${end_coord}" "${TARGET_DIR}${RUFUS_PRE_VCF}.gz" >> $COMBINED_PRE_VCF
        fi

        # Advance start coordinate
        ((start_coord=end_coord+1))
    done
    echo "Done with chr${curr_chr}"
done

bgzip $COMBINED_VCF
bcftools index -t "${COMBINED_VCF}.gz"

bgzip $COMBINED_PRE_VCF
bcftools index -t "${COMBINED_PRE_VCF}.gz"
                                           