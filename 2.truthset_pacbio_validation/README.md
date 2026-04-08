# 2. Truth Set Validation with PacBio Long Reads

This directory contains the scripts and notebooks used for the **PacBio long-read validation** step of the COLO829BLT50 reference mutation set, as described in the manuscript *"Detecting low allele fraction mosaic mutations: A multi-institution benchmarking with cell line mixtures."*

## Overview

Although ultra-deep Illumina sequencing improves the quality of variant calls, platform-specific errors persist — especially in genomic regions where read alignment is imperfect or where copy-number alterations cause a rare germline variant to appear as somatic. To eliminate short-read-specific artifacts, each candidate COLO829-specific mutation from [`1.truthset_illumina`](../1.truthset_illumina) was evaluated for support in **PacBio HiFi long-read data** (~170× COLO829 and ~330× COLO829BL, ~500× combined).

A variant was considered **valid** if:

1. The alternate allele was supported by **at least two PacBio reads** in COLO829 (to avoid overly permissive single-read validation, since chance errors are expected at such combined depths), AND
2. The alternate allele was **absent in COLO829BL** PacBio data.

For indels, an additional filter was applied using COLO829BL Illumina data to remove variants lacking tumor specificity, because long-read sequencing can produce stutter-associated indel artifacts near homopolymers and simple repeats.

The final validated reference set comprises **44,005 SNVs, 961 insertions, and 1,098 deletions**.

## Files

- **`Snakefile`** — Snakemake workflow that orchestrates the PacBio pileup / read-support evaluation across candidate SNV and indel sites.
- **`pipe_indel_pu.sh`** — Shell pipeline for generating per-site PacBio pileups at candidate indel positions.
- **`annotate_indel_info-wSecondAlt.py`** — Python script that annotates each candidate indel with supporting-read counts, including a second alternate allele field to handle multi-allelic cases and stutter artifacts.
- **`Validate_snv.ipynb`** — Jupyter notebook applying the two-read support rule (alt present in COLO829, absent in COLO829BL) to classify SNV candidates as valid / invalid.
- **`Validate_indel.ipynb`** — Jupyter notebook performing analogous validation for indels, including the additional Illumina-based tumor-specificity filter to remove long-read stutter artifacts.

## Workflow Summary

1. Receive the unified candidate SNV and indel set from [`1.truthset_illumina`](../1.truthset_illumina).
2. Run the Snakemake pipeline to compute per-site PacBio pileups in both COLO829 and COLO829BL.
3. Annotate candidate indel sites with detailed allele-support information (`annotate_indel_info-wSecondAlt.py`).
4. Execute the `Validate_snv.ipynb` and `Validate_indel.ipynb` notebooks to apply the final validation criteria.
5. Output the validated mosaic SNV / indel reference set for use in downstream benchmarking and negative-control construction.

## Validation Results (summary)

- **SNVs**: median validation rate of ~92% (range 79–96%) across callers; most unvalidated candidates were rejected due to the absence of the alternate allele in PacBio COLO829 (~4.3%) or the presence of the alternate allele in PacBio COLO829BL (~2.8%).
- **Indels**: only ~63% of candidates validated (range 48–91%); most unsupported indels were rejected because the tumor long-read data did not support them.

## Input Data

- Illumina candidate VCFs from [`1.truthset_illumina`](../1.truthset_illumina)
- **COLO829** PacBio HiFi WGS (~170× coverage)
- **COLO829BL** PacBio HiFi WGS (~330× coverage)

Raw sequencing data are publicly available via the SMaHT Data Portal: <https://data.smaht.org>.

## Next Step

The validated SNVs and indels produced here form the **reference (positive-control) set** used in [`5.benchmarking`](../5.benchmarking), together with the negative-control set from [`3.negative_control`](../3.negative_control).
