
#!/bin/bash

input_vcf="/n/data1/hms/dbmi/park/julia/SMaHT/benchmarking/mpileup/pileupVAF/pileup_all_calls_for_pileupVCFcalc.norm.vcf"
output_tsv="/n/data1/hms/dbmi/park/julia/SMaHT/benchmarking/mpileup/pileupVAF/pileup_all_calls_for_pileupVCFcalc.norm.pileupVAF.tsv"

# Get sample names
samples=$(bcftools query -l $input_vcf | sed 's#.*/##; s/\.bam$//')

# Build header
header="CHROM\tPOS\tREF\tALT\tDP_INFO"
for s in $samples; do
    header="${header}\t${s}_DP\t${s}_ADref\t${s}_ADalt"
done
header="${header}\tDP_sum_samples\tADalt_sum\tVAF_INFO\tVAF_SAMPLES"

# Write header to output
echo -e "$header" > $output_tsv

# Extract: DP_INFO, DP and AD
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP[\t%DP\t%AD]\n' $input_vcf |
awk '{
    DP_INFO = $5
    total_dp_samples = 0
    total_alt = 0
    per_sample = ""

    # Loop over each sample: first DP, then AD
    for (i=6; i<=NF; i+=2) {
        dp = $i
        split($(i+1), counts, ",")
        refc = counts[1]
        altc = counts[2]
        total_dp_samples += dp
        total_alt += altc
        per_sample = per_sample sprintf("\t%s\t%s\t%s", dp, refc, altc)
    }

    vaf_info = (DP_INFO > 0) ? total_alt / DP_INFO : "NA"
    vaf_samples = (total_dp_samples > 0) ? total_alt / total_dp_samples : "NA"

    print $1, $2, $3, $4, DP_INFO per_sample, total_dp_samples, total_alt, vaf_info, vaf_samples
}' OFS="\t" >> $output_tsv

