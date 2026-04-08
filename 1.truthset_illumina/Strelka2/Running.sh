#!/bin/bash

DIR=/n/data1/hms/dbmi/park/jiny/SMaHT/COLO829/TruthSet
SCRIPT=$DIR/workflow
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
DataPath=$DIR/input/merge_three_bams/
IntervalPath=/n/data1/hms/dbmi/park/jiny/resources/WGS_Intervals_30m_bed
BAM_BL=$(ls $DataPath | egrep "bam|BL" |  egrep -v "bai|merge_three")
BAM_BL_LIST=(${BAM_BL// / })

for idx in ${!BAM_BL_LIST[@]};
do
	for i in {1..91} X Y; #INTERVALS
	do
		#bed compress and index
		#/n/app/htslib/1.14/bin/bgzip $IntervalPath/interval_$i".bed"
		#/n/app/htslib/1.14/bin/tabix $IntervalPath/interval_$i".bed.gz"

		Interval=$IntervalPath/interval_$i".bed.gz"
                Interval_idx=$i

		T_bam=$DIR/input/merge_three_bams/COLO829T_Ill_200X.bam
		BL_bam=${BAM_BL_LIST[idx]}
		IDX=`echo $BL_bam | awk -F"." '{print $1}'`
		echo $BL_bam
		echo $IDX
		echo $interval
		:<<END
		#Manta - config generation
		OutPath=$DIR/output/Manta/$IDX'_'$i
		LogPath=$DIR/log/Manta
		sbatch \
		-A park \
		-p medium \
		-t 80:00:00 \
		--mem 8G \
		-e $LogPath/MAT_$IDX'_'$i.err \
		-o $LogPath/MAT_$IDX'_'$i.log \
		$SCRIPT/pipe_Manta_config.sh \
		$DataPath $BL_bam $T_bam $REF $Interval $Interval_idx $OutPath $IDX $i
END
		#Strelka2
		OutPath=$DIR/output/Strelka2/$IDX'_'$i
		MantaPath=$DIR'/output/Manta/'$IDX'_'$i
		echo $MantaPath
                LogPath=$DIR/log/Strelka2
                sbatch \
                -A park \
                -p medium \
		-t 80:00:00 \
		--mem 8G \
                -e $LogPath/STK_$IDX'_'$i.err \
                -o $LogPath/STK_$IDX'_'$i.log \
                $SCRIPT/pipe_Strelka2.sh \
                $DataPath $BL_bam $T_bam $REF $Interval $OutPath $IDX $i $MantaPath


		sleep 1
	done	
	sleep 1
done
