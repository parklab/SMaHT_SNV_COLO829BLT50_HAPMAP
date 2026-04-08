#!/bin/bash

#!/bin/bash

mkdir -p /n/data1/hms/dbmi/park/jiny/resources/WGS_Intervals_5m
mkdir -p /n/data1/hms/dbmi/park/jiny/resources/WGS_Intervals_5m_bed
python /n/data1/hms/dbmi/park/jiny/resources/scripts/Make_intervals_chrom.py \
-d /n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.dict \
-o /n/data1/hms/dbmi/park/jiny/resources/WGS_Intervals_5m \
-m 5 \
-ob /n/data1/hms/dbmi/park/jiny/resources/WGS_Intervals_5m_bed


