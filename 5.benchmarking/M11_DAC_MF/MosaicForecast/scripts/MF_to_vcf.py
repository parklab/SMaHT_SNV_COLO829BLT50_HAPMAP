#!/bin/python

import sys

# Define the input file path based on the provided TARGET_DIR and DIR_NAME

INPUT_FILE = sys.argv[1]
OUTPUT_FILE = sys.argv[2] #f"{INPUT_FILE}.vcf"

# VCF header
vcf_header = """##fileformat=VCFv4.2
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">
##contig=<ID=chr1,length=248956422,assembly=hg38>
##contig=<ID=chr2,length=242193529,assembly=hg38>
##contig=<ID=chr3,length=198295559,assembly=hg38>
##contig=<ID=chr4,length=190214555,assembly=hg38>
##contig=<ID=chr5,length=181538259,assembly=hg38>
##contig=<ID=chr6,length=170805979,assembly=hg38>
##contig=<ID=chr7,length=159345973,assembly=hg38>
##contig=<ID=chr8,length=145138636,assembly=hg38>
##contig=<ID=chr9,length=138394717,assembly=hg38>
##contig=<ID=chr10,length=133797422,assembly=hg38>
##contig=<ID=chr11,length=135086622,assembly=hg38>
##contig=<ID=chr12,length=133275309,assembly=hg38>
##contig=<ID=chr13,length=114364328,assembly=hg38>
##contig=<ID=chr14,length=107043718,assembly=hg38>
##contig=<ID=chr15,length=101991189,assembly=hg38>
##contig=<ID=chr16,length=90338345,assembly=hg38>
##contig=<ID=chr17,length=83257441,assembly=hg38>
##contig=<ID=chr18,length=80373285,assembly=hg38>
##contig=<ID=chr19,length=58617616,assembly=hg38>
##contig=<ID=chr20,length=64444167,assembly=hg38>
##contig=<ID=chr21,length=46709983,assembly=hg38>
##contig=<ID=chr22,length=50818468,assembly=hg38>
##contig=<ID=chrX,length=156040895,assembly=hg38>
##contig=<ID=chrY,length=57227415,assembly=hg38>
##contig=<ID=chrM,length=16569,assembly=hg38>
##contig=<ID=chr1_KI270706v1_random,length=175055,assembly=hg38>
##contig=<ID=chr1_KI270707v1_random,length=32032,assembly=hg38>
##contig=<ID=chr1_KI270708v1_random,length=127682,assembly=hg38>
##contig=<ID=chr1_KI270709v1_random,length=66860,assembly=hg38>
##contig=<ID=chr1_KI270710v1_random,length=40176,assembly=hg38>
##contig=<ID=chr1_KI270711v1_random,length=42210,assembly=hg38>
##contig=<ID=chr1_KI270712v1_random,length=176043,assembly=hg38>
##contig=<ID=chr1_KI270713v1_random,length=40745,assembly=hg38>
##contig=<ID=chr1_KI270714v1_random,length=41717,assembly=hg38>
##contig=<ID=chr2_KI270715v1_random,length=161471,assembly=hg38>
##contig=<ID=chr2_KI270716v1_random,length=153799,assembly=hg38>
##contig=<ID=chr3_GL000221v1_random,length=155397,assembly=hg38>
##contig=<ID=chr4_GL000008v2_random,length=209709,assembly=hg38>
##contig=<ID=chr5_GL000208v1_random,length=92689,assembly=hg38>
##contig=<ID=chr9_KI270717v1_random,length=40062,assembly=hg38>
##contig=<ID=chr9_KI270718v1_random,length=38054,assembly=hg38>
##contig=<ID=chr9_KI270719v1_random,length=176845,assembly=hg38>
##contig=<ID=chr9_KI270720v1_random,length=39050,assembly=hg38>
##contig=<ID=chr11_KI270721v1_random,length=100316,assembly=hg38>
##contig=<ID=chr14_GL000009v2_random,length=201709,assembly=hg38>
##contig=<ID=chr14_GL000225v1_random,length=211173,assembly=hg38>
##contig=<ID=chr14_KI270722v1_random,length=194050,assembly=hg38>
##contig=<ID=chr14_GL000194v1_random,length=191469,assembly=hg38>
##contig=<ID=chr14_KI270723v1_random,length=38115,assembly=hg38>
##contig=<ID=chr14_KI270724v1_random,length=39555,assembly=hg38>
##contig=<ID=chr14_KI270725v1_random,length=172810,assembly=hg38>
##contig=<ID=chr14_KI270726v1_random,length=43739,assembly=hg38>
##contig=<ID=chr15_KI270727v1_random,length=448248,assembly=hg38>
##contig=<ID=chr16_KI270728v1_random,length=1872759,assembly=hg38>
##contig=<ID=chr17_GL000205v2_random,length=185591,assembly=hg38>
##contig=<ID=chr17_KI270729v1_random,length=280839,assembly=hg38>
##contig=<ID=chr17_KI270730v1_random,length=112551,assembly=hg38>
##contig=<ID=chr22_KI270731v1_random,length=150754,assembly=hg38>
##contig=<ID=chr22_KI270732v1_random,length=41543,assembly=hg38>
##contig=<ID=chr22_KI270733v1_random,length=179772,assembly=hg38>
##contig=<ID=chr22_KI270734v1_random,length=165050,assembly=hg38>
##contig=<ID=chr22_KI270735v1_random,length=42811,assembly=hg38>
##contig=<ID=chr22_KI270736v1_random,length=181920,assembly=hg38>
##contig=<ID=chr22_KI270737v1_random,length=103838,assembly=hg38>
##contig=<ID=chr22_KI270738v1_random,length=99375,assembly=hg38>
##contig=<ID=chr22_KI270739v1_random,length=73985,assembly=hg38>
##contig=<ID=chrY_KI270740v1_random,length=37240,assembly=hg38>
##contig=<ID=chrUn_KI270302v1,length=2274,assembly=hg38>
##contig=<ID=chrUn_KI270304v1,length=2165,assembly=hg38>
##contig=<ID=chrUn_KI270303v1,length=1942,assembly=hg38>
##contig=<ID=chrUn_KI270305v1,length=1472,assembly=hg38>
##contig=<ID=chrUn_KI270322v1,length=21476,assembly=hg38>
##contig=<ID=chrUn_KI270320v1,length=4416,assembly=hg38>
##contig=<ID=chrUn_KI270310v1,length=1201,assembly=hg38>
##contig=<ID=chrUn_KI270316v1,length=1444,assembly=hg38>
##contig=<ID=chrUn_KI270315v1,length=2276,assembly=hg38>
##contig=<ID=chrUn_KI270312v1,length=998,assembly=hg38>
##contig=<ID=chrUn_KI270311v1,length=12399,assembly=hg38>
##contig=<ID=chrUn_KI270317v1,length=37690,assembly=hg38>
##contig=<ID=chrUn_KI270412v1,length=1179,assembly=hg38>
##contig=<ID=chrUn_KI270411v1,length=2646,assembly=hg38>
##contig=<ID=chrUn_KI270414v1,length=2489,assembly=hg38>
##contig=<ID=chrUn_KI270419v1,length=1029,assembly=hg38>
##contig=<ID=chrUn_KI270418v1,length=2145,assembly=hg38>
##contig=<ID=chrUn_KI270420v1,length=2321,assembly=hg38>
##contig=<ID=chrUn_KI270424v1,length=2140,assembly=hg38>
##contig=<ID=chrUn_KI270417v1,length=2043,assembly=hg38>
##contig=<ID=chrUn_KI270422v1,length=1445,assembly=hg38>
##contig=<ID=chrUn_KI270423v1,length=981,assembly=hg38>
##contig=<ID=chrUn_KI270425v1,length=1884,assembly=hg38>
##contig=<ID=chrUn_KI270429v1,length=1361,assembly=hg38>
##contig=<ID=chrUn_KI270442v1,length=392061,assembly=hg38>
##contig=<ID=chrUn_KI270466v1,length=1233,assembly=hg38>
##contig=<ID=chrUn_KI270465v1,length=1774,assembly=hg38>
##contig=<ID=chrUn_KI270467v1,length=3920,assembly=hg38>
##contig=<ID=chrUn_KI270435v1,length=92983,assembly=hg38>
##contig=<ID=chrUn_KI270438v1,length=112505,assembly=hg38>
##contig=<ID=chrUn_KI270468v1,length=4055,assembly=hg38>
##contig=<ID=chrUn_KI270510v1,length=2415,assembly=hg38>
##contig=<ID=chrUn_KI270509v1,length=2318,assembly=hg38>
##contig=<ID=chrUn_KI270518v1,length=2186,assembly=hg38>
##contig=<ID=chrUn_KI270508v1,length=1951,assembly=hg38>
##contig=<ID=chrUn_KI270516v1,length=1300,assembly=hg38>
##contig=<ID=chrUn_KI270512v1,length=22689,assembly=hg38>
##contig=<ID=chrUn_KI270519v1,length=138126,assembly=hg38>
##contig=<ID=chrUn_KI270522v1,length=5674,assembly=hg38>
##contig=<ID=chrUn_KI270511v1,length=8127,assembly=hg38>
##contig=<ID=chrUn_KI270515v1,length=6361,assembly=hg38>
##contig=<ID=chrUn_KI270507v1,length=5353,assembly=hg38>
##contig=<ID=chrUn_KI270517v1,length=3253,assembly=hg38>
##contig=<ID=chrUn_KI270529v1,length=1899,assembly=hg38>
##contig=<ID=chrUn_KI270528v1,length=2983,assembly=hg38>
##contig=<ID=chrUn_KI270530v1,length=2168,assembly=hg38>
##contig=<ID=chrUn_KI270539v1,length=993,assembly=hg38>
##contig=<ID=chrUn_KI270538v1,length=91309,assembly=hg38>
##contig=<ID=chrUn_KI270544v1,length=1202,assembly=hg38>
##contig=<ID=chrUn_KI270548v1,length=1599,assembly=hg38>
##contig=<ID=chrUn_KI270583v1,length=1400,assembly=hg38>
##contig=<ID=chrUn_KI270587v1,length=2969,assembly=hg38>
##contig=<ID=chrUn_KI270580v1,length=1553,assembly=hg38>
##contig=<ID=chrUn_KI270581v1,length=7046,assembly=hg38>
##contig=<ID=chrUn_KI270579v1,length=31033,assembly=hg38>
##contig=<ID=chrUn_KI270589v1,length=44474,assembly=hg38>
##contig=<ID=chrUn_KI270590v1,length=4685,assembly=hg38>
##contig=<ID=chrUn_KI270584v1,length=4513,assembly=hg38>
##contig=<ID=chrUn_KI270582v1,length=6504,assembly=hg38>
##contig=<ID=chrUn_KI270588v1,length=6158,assembly=hg38>
##contig=<ID=chrUn_KI270593v1,length=3041,assembly=hg38>
##contig=<ID=chrUn_KI270591v1,length=5796,assembly=hg38>
##contig=<ID=chrUn_KI270330v1,length=1652,assembly=hg38>
##contig=<ID=chrUn_KI270329v1,length=1040,assembly=hg38>
##contig=<ID=chrUn_KI270334v1,length=1368,assembly=hg38>
##contig=<ID=chrUn_KI270333v1,length=2699,assembly=hg38>
##contig=<ID=chrUn_KI270335v1,length=1048,assembly=hg38>
##contig=<ID=chrUn_KI270338v1,length=1428,assembly=hg38>
##contig=<ID=chrUn_KI270340v1,length=1428,assembly=hg38>
##contig=<ID=chrUn_KI270336v1,length=1026,assembly=hg38>
##contig=<ID=chrUn_KI270337v1,length=1121,assembly=hg38>
##contig=<ID=chrUn_KI270363v1,length=1803,assembly=hg38>
##contig=<ID=chrUn_KI270364v1,length=2855,assembly=hg38>
##contig=<ID=chrUn_KI270362v1,length=3530,assembly=hg38>
##contig=<ID=chrUn_KI270366v1,length=8320,assembly=hg38>
##contig=<ID=chrUn_KI270378v1,length=1048,assembly=hg38>
##contig=<ID=chrUn_KI270379v1,length=1045,assembly=hg38>
##contig=<ID=chrUn_KI270389v1,length=1298,assembly=hg38>
##contig=<ID=chrUn_KI270390v1,length=2387,assembly=hg38>
##contig=<ID=chrUn_KI270387v1,length=1537,assembly=hg38>
##contig=<ID=chrUn_KI270395v1,length=1143,assembly=hg38>
##contig=<ID=chrUn_KI270396v1,length=1880,assembly=hg38>
##contig=<ID=chrUn_KI270388v1,length=1216,assembly=hg38>
##contig=<ID=chrUn_KI270394v1,length=970,assembly=hg38>
##contig=<ID=chrUn_KI270386v1,length=1788,assembly=hg38>
##contig=<ID=chrUn_KI270391v1,length=1484,assembly=hg38>
##contig=<ID=chrUn_KI270383v1,length=1750,assembly=hg38>
##contig=<ID=chrUn_KI270393v1,length=1308,assembly=hg38>
##contig=<ID=chrUn_KI270384v1,length=1658,assembly=hg38>
##contig=<ID=chrUn_KI270392v1,length=971,assembly=hg38>
##contig=<ID=chrUn_KI270381v1,length=1930,assembly=hg38>
##contig=<ID=chrUn_KI270385v1,length=990,assembly=hg38>
##contig=<ID=chrUn_KI270382v1,length=4215,assembly=hg38>
##contig=<ID=chrUn_KI270376v1,length=1136,assembly=hg38>
##contig=<ID=chrUn_KI270374v1,length=2656,assembly=hg38>
##contig=<ID=chrUn_KI270372v1,length=1650,assembly=hg38>
##contig=<ID=chrUn_KI270373v1,length=1451,assembly=hg38>
##contig=<ID=chrUn_KI270375v1,length=2378,assembly=hg38>
##contig=<ID=chrUn_KI270371v1,length=2805,assembly=hg38>
##contig=<ID=chrUn_KI270448v1,length=7992,assembly=hg38>
##contig=<ID=chrUn_KI270521v1,length=7642,assembly=hg38>
##contig=<ID=chrUn_GL000195v1,length=182896,assembly=hg38>
##contig=<ID=chrUn_GL000219v1,length=179198,assembly=hg38>
##contig=<ID=chrUn_GL000220v1,length=161802,assembly=hg38>
##contig=<ID=chrUn_GL000224v1,length=179693,assembly=hg38>
##contig=<ID=chrUn_KI270741v1,length=157432,assembly=hg38>
##contig=<ID=chrUn_GL000226v1,length=15008,assembly=hg38>
##contig=<ID=chrUn_GL000213v1,length=164239,assembly=hg38>
##contig=<ID=chrUn_KI270743v1,length=210658,assembly=hg38>
##contig=<ID=chrUn_KI270744v1,length=168472,assembly=hg38>
##contig=<ID=chrUn_KI270745v1,length=41891,assembly=hg38>
##contig=<ID=chrUn_KI270746v1,length=66486,assembly=hg38>
##contig=<ID=chrUn_KI270747v1,length=198735,assembly=hg38>
##contig=<ID=chrUn_KI270748v1,length=93321,assembly=hg38>
##contig=<ID=chrUn_KI270749v1,length=158759,assembly=hg38>
##contig=<ID=chrUn_KI270750v1,length=148850,assembly=hg38>
##contig=<ID=chrUn_KI270751v1,length=150742,assembly=hg38>
##contig=<ID=chrUn_KI270752v1,length=27745,assembly=hg38>
##contig=<ID=chrUn_KI270753v1,length=62944,assembly=hg38>
##contig=<ID=chrUn_KI270754v1,length=40191,assembly=hg38>
##contig=<ID=chrUn_KI270755v1,length=36723,assembly=hg38>
##contig=<ID=chrUn_KI270756v1,length=79590,assembly=hg38>
##contig=<ID=chrUn_KI270757v1,length=71251,assembly=hg38>
##contig=<ID=chrUn_GL000214v1,length=137718,assembly=hg38>
##contig=<ID=chrUn_KI270742v1,length=186739,assembly=hg38>
##contig=<ID=chrUn_GL000216v2,length=176608,assembly=hg38>
##contig=<ID=chrUn_GL000218v1,length=161147,assembly=hg38>
##contig=<ID=chrEBV,length=171823,assembly=hg38>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""
chromosomes = [ f"chr{i}" for i in list(range(1,23))] + ['X', 'Y']
vcf_records = []

# Read the input file and process each line
with open(INPUT_FILE, 'r') as file:
    for line in file:
        parts = line.strip().split()
        if parts[0] != 'id':
            info_parts = parts[0].split("~")
            chrom = info_parts[1]
            pos = info_parts[2]
            ref = info_parts[3]
            alt = info_parts[4]
            qual = parts[1]
            af = parts[23]
            info = f"TYPE={parts[3]}"
            format = "AF"
            sample = af                
            vcf_record = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\tPASS\t{info}\t{format}\t{sample}"
            vcf_records.append(vcf_record)

vcf_output = vcf_header + "\n".join(vcf_records)

# Write the VCF output to a file
with open(OUTPUT_FILE + "_tmp", "w") as vcf_file:
    vcf_file.write(vcf_output)

print(f"VCF file has been created: {OUTPUT_FILE+ "_tmp"}")

import subprocess
subprocess.run(["/n/data1/hms/dbmi/park/dominika/testing/smaht/github/bioinfo/tools/sort-bgzip-vcf.sh", OUTPUT_FILE + "_tmp", OUTPUT_FILE])
