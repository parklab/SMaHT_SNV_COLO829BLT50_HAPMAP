import sys
from granite.lib import vcf_parser


# This script compares pileups between PacBio COLO829BL/COLO829T (q=1 and Q=1) and PacBio COLO829BL/COLO829T (q=30 and Q=30). Positions that appear to have reported an alternate allele  in q/Q=1 but they are not reported after filtering lower qualities and have no alleternate allele in q/Q=30 are sent further to examine in Illumina data to cross-check.
def load_variants(vcf):
    vcf_obj = vcf_parser.Vcf(vcf)
    variants = {}
    for vnt in vcf_obj.parse_variants():
        if f"{vnt.CHROM}_{vnt.POS}_{vnt.REF}" in variants.keys():
            variants[f"{vnt.CHROM}_{vnt.POS}_{vnt.REF}"].append(vnt.ALT)
        else:
            variants[f"{vnt.CHROM}_{vnt.POS}_{vnt.REF}"] = [vnt.ALT]
    return variants


vcf = sys.argv[1] #Mpilup q/Q=30 pileups VCF
index = sys.argv[2] #Job ID

vcf_obj = vcf_parser.Vcf(vcf)
genotype1 = "/n/data1/hms/dbmi/park-smaht_dac/analysis/GCC_UW/COLO829T/pacbioHifi_bulkFiberSeq/seq_data/COLO829T_Hifi_ReplacedRG.bam"
genotype2 = "/n/data1/hms/dbmi/park/dominika/testing/smaht/negative_control/WGS/PacBio/COLO829BL/multiple_rg_2/COLO829_ReplacedreadGroups.bam"
PacBioq1_variants = load_variants(
    f"PU1_PacBio/PU1_PacBio_{index}.norm.vcf"
)
alternates = {}

for vnt in vcf_obj.parse_variants():
    if vnt.ALT != "<*>":
        if f"{vnt.CHROM}_{vnt.POS}" not in alternates.keys():
            alternates[f"{vnt.CHROM}_{vnt.POS}"] = [vnt.ALT]
        else:
            alternates[f"{vnt.CHROM}_{vnt.POS}"].append(vnt.ALT)

with open(
    f"NegativeControlSet/intermediate/Step2_PacBio_PU30_ToCheckinIllumina/Step2_ToCheckInIllumina_{index}.bed",
    "w",
) as fo:
    for vnt in vcf_obj.parse_variants():
        if vnt.ALT == "<*>" and f"{vnt.CHROM}_{vnt.POS}" not in alternates.keys():
            AD_1 = int(vnt.get_genotype_value(genotype1, "AD").split(",")[-1])
            AD_2 = int(vnt.get_genotype_value(genotype2, "AD").split(",")[-1])

            if AD_1 == 0 or AD_2 == 0:
                for alt in PacBioq1_variants[f"{vnt.CHROM}_{vnt.POS}_{vnt.REF}"]:
                    if alt != "<*>":
                        fo.write(f"{vnt.CHROM}\t{int(vnt.POS)}\t{vnt.REF}\t{alt}\n")
