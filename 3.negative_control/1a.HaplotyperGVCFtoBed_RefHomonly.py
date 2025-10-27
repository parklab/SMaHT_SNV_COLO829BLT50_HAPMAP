from granite.lib import vcf_parser
import subprocess
import sys

vcf = vcf_parser.Vcf(sys.argv[1])
output_bed = sys.argv[2]
output_bed_tmp = output_bed + "tmp"
with open(output_bed_tmp, "w") as fo:
    for vnt in vcf.parse_variants():
        if vnt.ALT == "<NON_REF>":
            start = int(vnt.POS)
            end = int(vnt.get_tag_value("END"))
            fo.write(f"{vnt.CHROM}\t{start - 1}\t{end}\n")
            

subprocess.run(
    f"bedtools merge -i {output_bed_tmp} > {output_bed}", 
    shell=True, check=True
)


