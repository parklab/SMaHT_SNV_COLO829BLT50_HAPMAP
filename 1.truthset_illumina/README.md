# 1. Truth Set Generation from Illumina Short-Read Data

This directory contains scripts for the **initial candidate variant identification** step in constructing the high-confidence low-VAF somatic mutation reference set for the COLO829BLT50 cell line mixture, as described in the manuscript *"Detecting low allele fraction mosaic mutations: A multi-institution benchmarking with cell line mixtures."*

## Overview

To maximize sensitivity in initial discovery of tumor-specific mutations, candidate somatic SNVs and indels were called from ultra-deep Illumina short-read WGS of COLO829 (tumor, ~370×) and COLO829BL (matched normal, ~280×) using **four complementary variant callers** that span diverse detection principles. Variants identified by any of these methods were merged into a unified candidate set comprising **61,649 SNVs and 4,475 indels**, which was subsequently validated with PacBio long-read data (see [`2.truthset_pacbio_validation`](../2.truthset_pacbio_validation)).

## Variant Callers

Each subdirectory contains the pipeline / scripts used to run one of the four Illumina-based callers:

- **[`Mutect2/`](Mutect2)** — GATK Mutect2, a Bayesian somatic genotyper that uses local haplotype assembly. Implemented as a Snakemake workflow with supporting scripts.
- **[`Strelka2/`](Strelka2)** — Strelka2, which combines a probabilistic somatic model with adaptive indel error estimation and haplotype modeling. Run together with Manta for candidate indel detection (`pipe_Manta_config.sh`, `pipe_Strelka2.sh`, `Running.sh`).
- **[`VarNet/`](VarNet)** — A deep-learning-based image classifier for somatic variants. Includes duplicate removal (`1.RemoveDup.sh`) and the main VarNet calling script (`2.Varnet.sh`).
- **RUFUS** *(not included here)* — A k-mer-based detector that is free from reference bias. Unlike the other tools, RUFUS was run with a matched control; it was provided by a collaborating group and is therefore not re-implemented in this repository.

## Utility Scripts

- **`Make_intervals_chrom.py`** — Splits the genome into per-chromosome interval files used for parallelized variant calling.
- **`Make_intervals_chrom_run.sh`** — Wrapper shell script to execute the interval generation.

## Workflow Summary

1. Prepare per-chromosome intervals (`Make_intervals_chrom*`).
2. Run Mutect2, Strelka2, and VarNet on COLO829 (tumor) vs. COLO829BL (matched normal) Illumina data in single-sample mode (except RUFUS, which uses the matched control).
3. Retain only `PASS` calls from each caller.
4. Merge the results into a unified candidate set for downstream PacBio validation.

## Input Data

- **COLO829** Illumina WGS (~370× coverage)
- **COLO829BL** Illumina WGS (~280× coverage)

Raw sequencing data are publicly available via the SMaHT Data Portal: <https://data.smaht.org>.

## Next Step

The candidate variants produced here are validated against PacBio HiFi long-read data in [`2.truthset_pacbio_validation`](../2.truthset_pacbio_validation), and further filtered using the matched normal in [`3.negative_control`](../3.negative_control) to produce the final reference variant set.
