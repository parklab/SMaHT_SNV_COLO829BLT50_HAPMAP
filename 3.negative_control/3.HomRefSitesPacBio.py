import sys
from granite.lib import vcf_parser

vcf = sys.argv[1]
index = sys.argv[2]

vcf_obj = vcf_parser.Vcf(vcf)

genotype1 = "/n/data1/hms/dbmi/park-smaht_dac/analysis/GCC_UW/COLO829T/pacbioHifi_bulkFiberSeq/seq_data/COLO829T_Hifi_ReplacedRG.bam"
genotype2 = "/n/data1/hms/dbmi/park/dominika/testing/smaht/negative_control/WGS/PacBio/COLO829BL/multiple_rg_2/COLO829_ReplacedreadGroups.bam"
variants_output = set()
variant_remove = set()
with open(
    f"NegativeControlSet/Step1_{index}.bed",
    "w",
) as fo:
    with open(f"NegativeControlSet/intermediate/ForStep2_{index}.bed", "w",) as fo2:
        for vnt in vcf_obj.parse_variants():
            if vnt.ALT == "<*>":
                AD_1 = int(vnt.get_genotype_value(genotype1, "AD").split(",")[-1])
                AD_2 = int(vnt.get_genotype_value(genotype2, "AD").split(",")[-1])

                if AD_1 == 0 and AD_2 == 0:
                    variants_output.add(
                        f"{vnt.CHROM}\t{int(vnt.POS) - 1}\t{int(vnt.POS)}\n"
                    )
            else:
                variant_remove.add(f"{vnt.CHROM}\t{int(vnt.POS) - 1}\t{int(vnt.POS)}\n")

        for variant in variants_output:
            if variant not in variant_remove:
                fo.write(variant)
            else:
                fo2.write(variant)
