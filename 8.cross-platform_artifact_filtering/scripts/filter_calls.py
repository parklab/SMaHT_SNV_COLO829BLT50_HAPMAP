from granite.lib import vcf_parser
import argparse

def load_variants_control(vcf, control_bam):

    vcf_calls_control = vcf_parser.Vcf(vcf)
    variants = set()
    for vnt in vcf_calls_control.parse_variants():
        vnt_id = f"{vnt.CHROM}_{vnt.POS}_{vnt.REF}_{vnt.ALT}"

        AD = int(vnt.get_genotype_value(control_bam, "AD").split(",")[-1])
        if AD > 1:
            variants.add(vnt_id)
    return variants

def load_variants_target(vcf, target_bam):

    vcf_calls_control = vcf_parser.Vcf(vcf)
    variants = set()
    for vnt in vcf_calls_control.parse_variants():
        vnt_id = f"{vnt.CHROM}_{vnt.POS}_{vnt.REF}_{vnt.ALT}"

        AD = int(vnt.get_genotype_value(target_bam, "AD").split(",")[-1])
        if AD > 1:
            variants.add(vnt_id)
    return variants

parser = argparse.ArgumentParser(description="Filtartion of recurrent mutations")

parser.add_argument("--input_vcf_pu", help="PU VCF", required=True)
parser.add_argument("--target_bam", help="Target BAM", required=True)
parser.add_argument("--control_bam", help="Control BAM", required=True)
parser.add_argument("--input_vcf_calls", help="Calls VCF", required=True)
parser.add_argument("--output_vcf_calls", help="Output calls VCF", required=True)


args = parser.parse_args()

variants_pu_control = load_variants_control(args.input_vcf_pu, args.control_bam)
variants_pu_target = load_variants_target(args.input_vcf_pu, args.target_bam)
with open(args.output_vcf_calls, "w") as final_filtered_vcf: 
    vcf_calls = vcf_parser.Vcf(args.input_vcf_calls)

    
    vcf_calls.write_header(final_filtered_vcf)
    for vnt in vcf_calls.parse_variants():
        vnt_id = f"{vnt.CHROM}_{vnt.POS}_{vnt.REF}_{vnt.ALT}"
        if vnt_id in variants_pu_target and vnt_id not in variants_pu_control:
            vcf_calls.write_variant(final_filtered_vcf, vnt)