import sys
from granite.lib import vcf_parser

illumina_vcf = sys.argv[1]
index = sys.argv[2]

def load_variants(vcf):
    vcf_obj = vcf_parser.Vcf(vcf)
    variants = []
    for vnt in vcf_obj.parse_variants():
        if vnt.ALT in ["A", "C", "G", "T"]:
            variants.append(f"{vnt.CHROM}_{vnt.POS}_{vnt.REF}_{vnt.ALT}")
    return set(variants)

PacBio_q1_variants = []

with open(
    f"NegativeControlSet/intermediate/Step2_PacBio_PU30_ToCheckinIllumina/Step2_ToCheckInIllumina_{index}.bed",
    "r",
) as f:
    for line in f:
        PacBio_q1_variants.append(line.strip().replace("\t", "_"))

illumina_q1_variants = load_variants(
    f"{index}.vcf"
)

with open(
    f"NegativeControlSet/Step3_{index}.bed",
    "w",
) as fo:
    for variant in PacBio_q1_variants:
        if variant not in illumina_q1_variants:
            chrom, pos, _, _ = variant.split("_")
            fo.write(f"{chrom}\t{int(pos) -1}\t{pos}\n")
