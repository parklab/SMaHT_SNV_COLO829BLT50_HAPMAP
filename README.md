# SMaHT SNV and Indel Detection Benchmarking Study
This repository contains analysis scripts for the manuscript: "Detecting low allele fraction mosaic mutations: A multi-institution benchmarking with a cancer cell line mixture."
This work is part of the Somatic Mosaicism across Human Tissues (SMaHT) Network effort to establish best practices for identifying somatic variants with low allele frequency in healthy tissues. 

## Introduction
We designed a comprehensive benchmarking study to evaluate and address the challenges of detecting somatic mutations at low and ultra-low VAFs, focusing on variants with VAF below 3%, ranging from below 1% to over 15%.

### Experimental Design
We employed two experimental designs to artificially generate samples with low VAF variants with known ground-truth:
- COLO829BLT50 Mixture: A biologically controlled mixture of 2% COLO829 (melanoma cell line) and 98% COLO829BL (matched normal) to generate somatic variants below 3% VAF. 
- HapMap Mixture: A mixture of six well-characterized HapMap samples at various fractions to generate artificial somatic variants with VAFs ranging from as low as 0.25% up to 16.5%.

To ensure a high-confidence benchmark, we combined diverse sequencing technologies and specialized computational approaches to generate reliable reference variant sets, benchmark somatic variant calling approaches, and provide  recommendations.

## Overview
The benchmarking process included multiple stages, beginning with the generation of a reliable reference variant set, followed by comparison of somatic variant detection approaches.

Focusing on the COLO829BLT50 experiment, variants for the reference set were called with three variant callers from short read Illumina data of the tumor line COLO829, validated with long read PacBio data of the same line, and filtered against variants from the matching normal line COLO829BL as negative control. Scripts for these steps can be found in 
[1.truthset_illumina](1.truthset_illumina), 
[2.truthset_pacbio_validation](2.truthset_pacbio_validation),
[3.negative_control](3.negative_control), respectively.
The reference sets for 
[SNVs](Resource/Truthset_SNV_COLO829.txt)
 and 
[Indels](Resource/Truthset_Indel_COLO829.txt)
as well as the 
[culture variants](Resource/Truthset_CultureDerived_SNV_COLO829.txt)
 are provided in the [Resource](Resource) directory.

As the cell lines were cultured in preparation for the generation of the COLO829BLT50 mixture, they acquired *de novo* variants that are absent from the initial COLO829 tumor and COLO829BL matched normal cell line. How those admixture-specific variants were identified is described in 
[4.admixture-only](4.admixture-only).

We benchmarked 12 different somatic variant calling approaches against the thoroughly prepared reference set of variants. All approaches are described in detail in the accompanying publication. While 11 approaches were provided by other groups in the SMaHT consortia, we employed Mutect2 and MosaicForecast in our provided approach. Details are described in 
[5.benchmarking](5.benchmarking). Additionally, to mitigate recurrent, alignment-driven false positive SNVs enriched in difficult and extreme genomic regions, we implemented a cross-sample pileup filtering workflow using an independent control sample, which we recommend as best practice. Details for this recommended filtering step are described in 
[8.cross-platform_artifact_filtering](8.cross-platform_artifact_filtering). 
The resulting call sets of all 12 approaches can be found in [7.call_sets](7.call_sets).

The reference variant set for the HapMap mixture experiment was generated from publicly available data from the 6 well-characterized cell lines as described in [HapMap_reference_variants](6.hapmap/HapMap_reference_variants/) and provided as [VCF](Resource/HapMapMixture_variant_reference_set.vcf.gz). [Sampling](6.hapmap/HapMap_somatic_variant_sampling_and_calling/sampling/)
as well as
[variant calling approaches](6.hapmap/HapMap_somatic_variant_sampling_and_calling/calling/) 
for the HapMap mixture are described in [6.hapmap](6.hapmap). 

## Cross-platform artifact filtering 


The [cross-platform artifact filtering  workflow](8.cross-platform_artifact_filtering/Snakemake/Snakefile) implements a cross-sample artifact filtering approach for short-read SNV calls. Candidate variant sites are evaluated using targeted pileups generated with bcftools mpileup, and allele depth (AD) is compared between a target sample and an independent control sample. The pipeline is implemented in Snakemake and outputs a filtered VCF of high-confidence variants.
