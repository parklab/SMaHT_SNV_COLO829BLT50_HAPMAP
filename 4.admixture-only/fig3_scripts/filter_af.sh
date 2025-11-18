#!/bin/bash

VCF=$1
OUT_VCF=$2
AF_MIN=$3
AF_MAX=$4

bcftools filter -i "INFO/ALLELIC_PCT>=$AF_MIN && INFO/ALLELIC_PCT<$AF_MAX" $VCF -o $OUT_VCF -Ov
