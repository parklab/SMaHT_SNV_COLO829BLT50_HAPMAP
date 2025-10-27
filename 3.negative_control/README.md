# Negative Control for COLO829BLT Benchmark


To define the negative control set, we first identified homozygous reference sites shared by  COLO829 and COLO829BL Illumina WGS. To define potential non-variant positions shared by COLO829 and COLO829BL, Santieon Haplotyper30 in the GVCF mode was applied independently on the two samples and putative homozygous positions shared by both were intersected for further processing. Next, we applied bcftools mpileup on the alignment of COLO829 and COLO829BL PacBio WGS in the putative homozygous candidates with the minimum mapping quality for an alignment (q parameter) and minimum base quality for a base to be considered (Q parameter) set to 1. Positions where no alternate allele found in both COLO829 and COLO829BL were collected as homozygous reference sites and included in the negative control set. The rest of the positions were then subjected to the next iteration of the mpileup analysis but with increased q and Q parameters (both set to 30). If alternate allele wasobserved with the higher quality threshold, those were removed from the negative control. Then, we also checked Illumina WGS data for to confirm remaining ambiguity. We again utilizedbcftools mpileup on the short-reads with adjusted q and Q parameters (set to 1) and juxtaposedwith the results of PacBio WGS pileups. Variants supported by both platforms were excluded, whereas those lacking Illumina support were retained as high-confidence negatives.

## Scripts

`1.Haplotyper_GVCF.sh`: Sentieon Haplotyper (gatk HaplotypeCaller equivalent) in GVCF mode 

`2.Mpileup.sh`: Bcftools mpileup for 2 BAMs COLO829BL and COLO829T with parametrized `--min-bq` and `--min-mq`

`3.HomRefSitesPacBio.py`: Script for parsing the output VCF file from `bcftools mpileup` on COLO829T and COLO829BL PacBio, positions with no alletrate allele AD=0 in both COLO829T and COLO829BL are incorporated into the Negative Control. Positions with alteranate allele in either sample are saved to a BED file for further evaluation.

`4.IlluminaCrossCheckSites.py`: Script for evaluation of PacBio COLO829T and COLO829BL pileups with `--min-bq=1` and `--min-mq=1` vs. --min-bq=30 and --min-mq=30. Sites with reported variants in both are removed from the Negative Control. Variants reported in only in the pileups the minimim qualities of 1, but are not reported in higher mapping abse qualities are passed for further evaluation. 

`5.CrossCheckPacBioIllumina.py`: Sites reported to be cross-checked between PacBio and Illumina are evaluated whether they occur in any Illumina sample. 