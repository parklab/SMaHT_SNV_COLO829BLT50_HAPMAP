# SNV Benchmark analysis on HapMap cell line mixture 
SMaHT designed the HapMap Mixture Benchmarking experiment to evaluate somatic variant calling across a large range of VAFs. Six well-characterized HapMap cell lines (HG002, HG005, HG00438, HG02257, HG02486, and HG02622) were mixed to create target VAFs ranging from 0.25% to 16.5%. 

The SNV analysis of this benchmarking experiment comprised:
1. the generation of a reference set of variants 
2. the detection (sampling) and calling of somatic variants in sequencing data of the cell line mixture

## Generation of HapMap Mixture reference SNV set
Publicly available variants as well as confident regions of the 6 HapMap cell lines were obtained from GIAB and Pangenome Projects, 
as shown in [get_and_filter_variants.sh](HapMap_reference_variants/get_and_filter_variants.sh) and [get_confident_regions.sh](HapMap_reference_variants/get_confident_regions.sh). 
Variants were filtered, combined and expected VAFs were calculated according to the fractions of cell lines in the mixture, using 
[HapMap_reference_set_expected_VAF.R](HapMap_reference_variants/HapMap_reference_set_expected_VAF.R). 

## SNV detection

### Subsampling of sequencing data
The HapMap mixture was sequenced across multiple centers at varying target depths, reaching up to 470x, to evaluate the impact of coverage on variant detection sensitivity and accuracy.
High-coverage data was subsampled to generate multiple 100x datasets using 
[downsample_bam.sh](HapMap_somatic_variant_sampling_and_calling/subset_and_merge_samples/downsample_bam.sh). 
To evaluate if variant detection differs between a single high-coverage dataset and multiple lower-coverage datasets merged together, we generated merged files using 
[merge_bam.sh](HapMap_somatic_variant_sampling_and_calling/subset_and_merge_samples/merge_bam.sh). 

### Sampling of reference variants
We first checked whether refence variants were present in the sequencing data of the mixture by creating pileups on target postions 
[1.run_pileup_on_target_positions.sh](HapMap_somatic_variant_sampling_and_calling/sampling/1.run_pileup_on_target_positions.sh) in genomic intervals, which were merged and normalized using 
[2.merge_and_norm_pileup.sh](HapMap_somatic_variant_sampling_and_calling/sampling/2.merge_and_norm_pileup.sh). 
In order to compare to the expected VAFs calculated from the mixture ratio of the cell lines, the observed VAFs were calculated from the pileup using 
[3.calculate_pileupVAF.sh](HapMap_somatic_variant_sampling_and_calling/sampling/3.calculate_pileupVAF.sh).

### Calling of somatic variants
#### short read Illumina sequencing data
Somatic variants were called from short read data with TNHaplotyper2 by running 
[1.run_TNhaplotyper2.sh](HapMap_somatic_variant_sampling_and_calling/calling/Illumina_TNHaplotyper/1.run_TNhaplotyper2.sh)
and
[2.run_TNFilter.sh](HapMap_somatic_variant_sampling_and_calling/calling/Illumina_TNHaplotyper/2.run_TNFilter.sh).
The mixture of 6 cell lines from 6 different individuals resulted in a significantly higher number of target variants than typically expected in normal samples.
Additionally, these variants are distributed across 12 distinct haplotypes, which is highly untypical for normal samples.
To acknowledge this artifical design, variants flagged as “clustered_events” (indicating multiple events on the same haplotype)
or “haplotype” (variants co-occurring with filtered variants on the same haplotype) were retained alongside PASS calls.

#### long read PacBio sequencing data 
Somatic variants were called from long read data with ClairS-TO by running
[somatic_variant_calling_on_PacBio_w_ClairS_TO.sh](HapMap_somatic_variant_sampling_and_calling/calling/PacBio_ClairS_TO/somatic_variant_calling_on_PacBio_w_ClairS_TO.sh).
Variants were not filtered against a panel of normals to avoid removing true positives. 
Calls flagged as “NonSomatic” (indicating overlap with germline SNV databases such as 1000 Genomes, gnomAD, or dbSNP) were retained, 
as the HapMap donor germline variants are expected to appear in these databases.
