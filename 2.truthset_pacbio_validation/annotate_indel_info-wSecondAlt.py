#!/home/yoh855/miniconda3/bin/python
# -*-coding:utf-8 -*-
'''
@Time: 2024/05/17
@Author: Yu Wei Zhang
@Contact: yuwei_zhang@hms.harvard.edu
Jiny modified ispass and returning dic for q_seq
'''
from pysam import VariantFile, tabix_index, AlignmentFile
import os
from collections import defaultdict 
import re


class AnnSV:
    def __init__(self, in_vcf=None, 
                 bam_file=None, 
                 out_vcf=None,
                 min_mapq=1,
                 dist_threshold=51,
                 verbose=False):
        self.in_vcf = in_vcf
        self.bam_file = bam_file
        self.bam = None
        self.out_vcf = out_vcf
        self.dist_threshold = dist_threshold
        self.min_mapq = min_mapq
        self.verbose = verbose

    def _checkFile(self):
        if not os.path.exists(self.in_vcf):
            raise FileNotFoundError("Input vcf file not found.")

        if not os.path.exists(self.bam_file):
            raise FileNotFoundError("Bam file not found.")
        
        # Extract the file name of self.bam_file
        bam_file_name = os.path.basename(self.bam_file).replace(".bam", "")
        
        if self.out_vcf is None:
            self.out_vcf = self.in_vcf.replace(".vcf", "_annoAllAlt_{}.vcf".format(bam_file_name))
        print("Output vcf file: {}".format(self.out_vcf))

    def _ispass(self, record):
        if ('PASS' in list(record.filter) and list(record.filter) != []) or list(record.filter) == []:
            return True
        else:
            return True

    def run(self):
        self._checkFile()
        vcf = VariantFile(self.in_vcf)
        if not "AF" in vcf.header.info:
            vcf.header.info.add("AF", number=1, type="Float", description="Allele frequency of the variant")
        if not "RNAMES" in vcf.header.info:
            vcf.header.info.add("RNAMES", number=1, type="String", description="Read names supporting the variant")
        if not "Depth" in vcf.header.info:
            vcf.header.info.add("Depth", number=1, type="Integer", description="Number of reads at the breakpoint")
        if not "DR" in vcf.header.info:
            vcf.header.info.add("DR", number=1, type="Integer", description="Number of reads supporting reference at the breakpoint")
        if not "DTV" in vcf.header.info:
            vcf.header.info.add("DTV", number=1, type="Integer", description="Number of reads supporting the targeted alternative allele at the breakpoint")
        if not "OtherAlt" in vcf.header.info:
            vcf.header.info.add("OtherAlt", number=1, type="String", description="Number of reads supporting all types of alternative allele at the breakpoint")
        annoVcf = VariantFile(self.out_vcf, "w", header=vcf.header)
        self.bam = AlignmentFile(self.bam_file)
        for record in vcf.fetch():
            if self._ispass(record) and record.chrom in self.bam.references:
                supp_reads, ref_reads, cov, all_alts = self.get_supp_reads(record)
                if self.verbose:
                    print(record, cov, len(supp_reads), len(ref_reads))
                record.info["AF"] = len(supp_reads) / cov if cov != 0 else 0
                record.info["RNAMES"] = ",".join(supp_reads)
                record.info["Depth"] = cov
                record.info["DR"] = len(ref_reads)
                record.info["DTV"] = len(supp_reads)
                record.info["OtherAlt"] = all_alts
                annoVcf.write(record)
        self.bam.close()
        annoVcf.close()
        tabix_index(self.out_vcf, preset="vcf", force=True)

    def get_supp_reads(self, record):
        if self.verbose:
            print(record)
        chrom = record.chrom
        ref = record.ref
        alts = record.alts
        if len(alts) > 1:
            raise ValueError("Compound Het detected")
        alt = alts[0]
        r_pos = record.pos
        r_stop = record.pos + len(ref)
        supp_reads = set()  
        ref_reads = set()
        other_reads = set()
        
        dic_q_seq = {}
        for read in self.bam.fetch(chrom, r_pos, r_stop):
            if read.mapping_quality >= self.min_mapq:
                ifSupport, q_seq = self._checkIndelInAln(read, r_pos, r_stop, ref, alt)
                if ifSupport == "Target":
                    supp_reads.add(read.query_name)
                elif ifSupport == "Ref":
                    ref_reads.add(read.query_name)
                elif ifSupport == "Other":
                    other_reads.add(read.query_name)
                if q_seq:
                    if q_seq not in dic_q_seq:
                        dic_q_seq[q_seq] = 1
                    else:
                        dic_q_seq[q_seq] += 1
        #print (dic_q_seq)
        all_alts = self._format_alt_sequences(dic_q_seq)
        return supp_reads, ref_reads, len(supp_reads) + len(ref_reads) + len(other_reads), all_alts
    
    def _checkIndelInAln(self, read, r_pos, r_stop, ref, alt):
        isSupp = "Ref"
        q_seq = None
        mapper = read.get_aligned_pairs(matches_only=True)  # (query_pos, ref_pos)
        q_pos = q_stop = None
        for pair in mapper:
            if pair[1] == r_pos - 1:
                q_pos = pair[0]
            if pair[1] == r_stop - 1:
                q_stop = pair[0]
        if q_pos is None and q_stop is None:
            isSupp = None # invalid read alignment
        elif q_pos is None:   # left soft-clipped
            if q_stop == mapper[0][0] and read.query_sequence[q_stop-len(alt):q_stop] == alt:
                isSupp = "Target"
            else:
                isSupp = "Other"
        elif q_stop is None:  # right soft-clipped
            if q_pos == mapper[-1][0] and read.query_sequence[q_pos:q_pos+len(alt)] == alt:
                isSupp = "Target"
            else:
                isSupp = "Other"
        else:
            q_seq = read.query_sequence[q_pos:q_stop]
            if q_stop - q_pos == len(alt):
                if len(alt) == 1:   # del
                    isSupp = "Target"
                else:  # ins, further check the inserted seq
                    if q_seq == alt:
                        isSupp = "Target"
                    elif q_seq == ref:
                        isSupp = "Ref"
                    else:
                        isSupp = "Other"
            elif q_seq == ref:
                isSupp = "Ref"
            else:
                isSupp = "Other"
        if self.verbose:
            print("{} ----- read: {}, q_pos: {}, q_stop: {}, q_seq: {}, ref: {}, alt: {}".format(
                isSupp, read.query_name, q_pos, q_stop, read.query_sequence[q_pos:q_stop], ref, alt))
        return isSupp, q_seq

    def _format_alt_sequences(self, dic_q_seq):
        dic_q_seq_sorted = sorted(dic_q_seq, key=dic_q_seq.get, reverse=True)
        l_q_seq = []
        for alt in dic_q_seq_sorted:
            l_q_seq.append(alt + '-' + str(dic_q_seq[alt]))
        return '_'.join(l_q_seq)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Annotate SVs from given dataset.")
    parser.add_argument("-i", "--input", help="Input indel vcf file.", required=True)
    parser.add_argument("-b", "--bam", help="Bam file.", required=True)
    parser.add_argument("-o", "--output", help="Output vcf file with AF and RNAMES info added.")
    parser.add_argument("-q", "--min_mapq", help="Minimum mapping quality.", default=1, type=int)
    parser.add_argument("-d", "--dist_threshold", help="Distance threshold for read alignment.", default=51, type=int)
    parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
    args = parser.parse_args()
    anno = AnnSV(args.input, args.bam, args.output, args.min_mapq, args.dist_threshold, args.verbose)
    anno.run()
