import os
import argparse

outdir = "/n/data1/hms/dbmi/park/jiny/resources/WGS_Intervals"
outdir_bed = "/n/data1/hms/dbmi/park/jiny/resources/WGS_Intervals_bed"
#f_fasta_dict = open("/n/data1/hms/dbmi/park/jiny/resources/WGS_Intervals_chr22/hg38_no_alt_chr22.dict", "r")
f_fasta_dict = open("/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.dict", "r")
global order_chr 
order_chr = [str(chrom) for chrom in list(range(1,23)) + ["X", "Y"]]
#order_chr = [str(chrom) for chrom in list(range(22,23))]
#state_chrM = True
length_chr = dict()
alt_contigs_by_chr = dict()

#Input

#dic_fh, nbases, output_dir, no_chr_prefix
#dic_fh= 
#nbases = 30 * 10 ** 6
#output_dir= "/home/yoh855/workflow/TruthSet_running/output/Intervals"
#no_chr_prefix = "--no-chr-prefix"

def make_interval_file(f_fasta_dict, nbases, output_dir, output_dir_bed):
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    global order_chr
    order_chr = [f"chr{chrom}" for chrom in order_chr]
    
    #Put length of chromosomes in dict
    for line in f_fasta_dict:
        if line.startswith("@SQ"):
            try:
                d = dict([entry.split(":",1) for entry in line.strip().split("\t")[1:]])
            except:
                print (line)
                raise
            if d["SN"] in order_chr:
                length_chr[d["SN"]] = int(d["LN"])
            elif d["SN"].endswith("_alt"):
                alt_contigs_by_chr[d["SN"].split("_")[0][d["SN"]]] = int(d("LN"))

    current_length = 0
    current_index = 0

    sequences = order_chr[:]
    order = iter(sequences)
    current_chrom = next(order)
    current_position = 0

    fh= open(f"{output_dir}/interval_{current_index}.list", "w")
    fh2= open(f"{output_dir_bed}/interval_{current_index}.bed", "w")
    length_remaining_in_current_intervalfile = nbases
    while True:
        remaining_in_chrom = length_chr[current_chrom] - current_position
        current_position +=1
        if remaining_in_chrom >= length_remaining_in_current_intervalfile:
            ending_position = current_position + length_remaining_in_current_intervalfile + 1
            fh.write(f"{current_chrom}:{current_position}-{ending_position}\n")
            fh.close()
            fh2.write(f"{current_chrom}\t{current_position-1}\t{ending_position}\n")
            fh2.close()
            if remaining_in_chrom == length_remaining_in_current_intervalfile: #if we can write one more interval with current chrom
                try:
                    current_chrom = next(order)
                    current_position = 0
                    remaining_in_chrom = length_chr[current_chrom]
                except StopIteration: #If it was ended in the last chrom
                    done = True
                    break
            else: #if more to write on current chrom, in next interval file
                current_position = ending_position
            current_index +=1
            fh= open(f"{output_dir}/interval_{current_index}.list", "w")
            fh2= open(f"{output_dir_bed}/interval_{current_index}.bed", "w")
            length_remaining_in_current_intervalfile = nbases
        else: #if remaining is smaller, write the rest and 
            ending_position = current_position + remaining_in_chrom -1
            fh.write(f"{current_chrom}:{current_position}-{ending_position}\n")
            fh2.write(f"{current_chrom}\t{current_position-1}\t{ending_position}\n")
            try:
                current_chrom = next(order)
                current_position = 0
                remaining_in_chrom = length_chr[current_chrom]
            except StopIteration:
                break
    if not fh.closed:
        fh.close()
    if not fh2.closed:
        fh2.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = __doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-d",
        "--dict",
        type=argparse.FileType("r"),
        help="the .dict file to read"
    )
    parser.add_argument(
        "-m",
        "--megabases",
        type=int,
        help="Number of megabases to be included in each interval"
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        default=os.getcwd(),
        help="Output directory for interval files"
    )
    parser.add_argument(
        "-ob",
        "--output_directory_bed",
        default=os.getcwd(),
        help="Output directory for interval bed files"
    )
    args = parser.parse_args()
    make_interval_file(
        args.dict,
        args.megabases * 10 ** 6,
        args.output_directory,
        args.output_directory_bed
    )


