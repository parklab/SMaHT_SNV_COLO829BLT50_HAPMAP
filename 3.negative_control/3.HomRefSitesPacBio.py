import sys
from granite.lib import vcf_parser

vcf = sys.argv[1]
index = sys.argv[2]
tumor_bam_genotype = sys.argv[3]
blood_bam_genotype = sys.argv[4]

vcf_obj = vcf_parser.Vcf(vcf)

variants_output = list()
variant_remove = set()
with open(
    f"NegativeControlSet/Step1_{index}.bed",
    "w",
) as fo:
    with open(f"NegativeControlSet/intermediate/ForStep2_{index}.bed", "w",) as fo2:
        for vnt in vcf_obj.parse_variants():
            if vnt.ALT == "<*>":
                AD_1 = int(vnt.get_genotype_value(tumor_bam_genotype, "AD").split(",")[-1])
                AD_2 = int(vnt.get_genotype_value(blood_bam_genotype, "AD").split(",")[-1])

                if AD_1 == 0 and AD_2 == 0:
                    variants_output.append(
                        f"{vnt.CHROM}\t{int(vnt.POS) - 1}\t{int(vnt.POS)}\n"
                    )
            else:
                variant_remove.add(f"{vnt.CHROM}\t{int(vnt.POS) - 1}\t{int(vnt.POS)}\n")

        for variant in variants_output:
            if variant not in variant_remove:
                fo.write(variant)
            else:
                fo2.write(variant)
