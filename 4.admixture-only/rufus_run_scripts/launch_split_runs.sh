#!/bin/bash
#SBATCH --time=3-00:00:00                                                                                                          
#SBATCH --account=marth-rw 
#SBATCH --partition=marth-shared-rw 
#SBATCH -o /scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/rufus_runs/blt50_1mb_bcm/slurm_out/out_%j.out
#SBATCH -e /scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/rufus_runs/blt50_1mb_bcm/slurm_out/err_%j.err
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=u0746015@utah.edu 

# This is an example script of how all regional RUFUS jobs are launched.

# This script is used to piecemeal launch multiple RUFUS runs from start to finish. 
# If the number of jobs exceeds the user limit (hardcoded to 1000), will wait to launch subsequent
# jobs until more resources are available.

# Files and dirs for runscript
SUBJECT=/scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/BLT50/bcm.400x.merged.SMHTCOLO829BLT50.sorted.bam
CONTROL=/scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/indel_realigned_bams/COLO829BL_Ill_230X.bam
OUT_SLURM_DIR=/scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/rufus_runs/blt50_1mb_bcm/slurm_out/
REFERENCE=/scratch/ucgd/lustre/work/marth/shared/references/human/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
RESOURCE_DIR=/scratch/ucgd/lustre/work/marth/resources/RUFUS/
OUT_DIR=/scratch/ucgd/lustre-labs/marth/scratch/u0746015/COLO829/rufus_runs/blt50_1mb_bcm/
UNID=u0746015

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

SLURM_SCRIPT="#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --account=marth-rw
#SBATCH --partition=marth-rw
#SBATCH -o ${OUT_SLURM_DIR}%j/out.out
#SBATCH -e ${OUT_SLURM_DIR}%j/error.err
module load samtools/1.16
module load bamtools/2.5.1
module load bedtools/2.28.0
module load htslib/1.16
module load gcc/10.2.0
module load RUFUS/latest"

n=${#CHRS[@]}
for (( i = 0; i < n; i++ ))
do
    curr_len="${CHR_LENGTHS[i]}"
    curr_chr="${CHRS[i]}"
    ((chunks=curr_len/CHUNK_SIZE+1))
    echo "number of chunks is $chunks"
 
    start_coord=1
    echo "Starting loop for index $i chr${curr_chr}"
    for (( c = 1; c <= chunks; c++ ))
    do
	    curr_jobs=($(squeue -u $UNID | wc -l))
        while [ $curr_jobs -gt 980 ]; do
		    sleep 60;
	        curr_jobs=($(squeue -u $UNID | wc -l))
	    done

        ((end_coord=start_coord+CHUNK_SIZE-1))
        if [ $end_coord -gt $curr_len ]; then
            echo "passed end coordinate, setting to $curr_len"
            end_coord=$curr_len
        fi
	
	    # Make region dir
        TARGET_DIR="chr${curr_chr}/chr${curr_chr}_${start_coord}_${end_coord}/"
        mkdir -p $TARGET_DIR
        
	    # Compose slurm script
	    SLURM_FILE="${TARGET_DIR}blt50_500x.sh"
	    echo -e "$SLURM_SCRIPT" > $SLURM_FILE 
        echo -e "cd ${OUT_DIR}${TARGET_DIR} \n" >> $SLURM_FILE
	    echo -e "bash \$RUFUS_ROOT/runRufus.sh -s $SUBJECT -c $CONTROL -r $REFERENCE -f ${RESOURCE_DIR}GRCh38_full_analysis_set_plus_decoy_hla.25.Jhash -m 5 -k 25 -t 40 -L -vs -R chr${curr_chr}:${start_coord}-${end_coord}" >> $SLURM_FILE 
        
        # Get run going
	    sbatch $SLURM_FILE

        # Advance start coordinate
        ((start_coord=end_coord+1))
    done
    echo "Done with chr${curr_chr}"
done
