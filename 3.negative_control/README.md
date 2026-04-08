# Negative Control for COLO829BLT Benchmark

## Overview

To define the negative control set, we first identified homozygous reference sites shared by COLO829 and COLO829BL Illumina whole-genome sequencing (WGS).

To determine potential non-variant positions shared by COLO829 and COLO829BL, **Sentieon Haplotyper** (GATK HaplotypeCaller equivalent) was applied independently to both samples in **GVCF mode**. Putative homozygous reference positions reported in both samples were intersected and used as candidate sites for further processing.

Next, we applied `bcftools mpileup` on the **PacBio WGS alignments** of COLO829 and COLO829BL at these candidate sites. The minimum mapping quality for alignments (`-q`) and the minimum base quality (`-Q`) were initially set to **1**. Positions where **no alternate allele** was observed in either COLO829 or COLO829BL were collected as homozygous reference sites and included in the **negative control set**.

The remaining positions were subjected to a second iteration of the mpileup analysis with increased quality thresholds (`-q 30` and `-Q 30`). If an alternate allele was observed at these higher quality thresholds, those positions were removed from the negative control.

Finally, we evaluated the remaining ambiguous sites using **Illumina WGS data**. We again ran `bcftools mpileup` on the Illumina short reads with `-q 1` and `-Q 1` and compared the results with the PacBio pileups. Positions supported by **both sequencing platforms** were excluded, whereas positions lacking Illumina support were retained as **high-confidence negative sites**.

---

# Scripts

### `1.Haplotyper_GVCF.sh`
Runs **Sentieon Haplotyper** (equivalent to GATK HaplotypeCaller) in **GVCF mode**.

### `1a.HaplotyperGVCFtoBed_RefHomonly.py`
Extracts **homozygous reference regions** reported by Haplotyper in GVCF mode and converts them to **BED format**.

### `2.Mpileup.sh`
Runs `bcftools mpileup` for **COLO829BL and COLO829T BAM files** with configurable parameters:

- `--min-bq` (minimum base quality)
- `--min-mq` (minimum mapping quality)

### `3.HomRefSitesPacBio.py`
Parses the VCF output from `bcftools mpileup` on **PacBio COLO829T and COLO829BL**.

- Positions with **no alternate allele (AD=0)** in both samples are added to the **negative control**.
- Positions with alternate alleles in either sample are written to a **BED file** for further evaluation.

### `4.IlluminaCrossCheckSites.py`
Evaluates PacBio pileups generated with two quality thresholds:

- `--min-bq=1`, `--min-mq=1`
- `--min-bq=30`, `--min-mq=30`

Processing logic:

- Sites with variants detected under **both thresholds** are removed from the negative control.
- Sites detected only under **low-quality thresholds** are passed for additional evaluation.

### `5.CrossCheckPacBioIllumina.py`
Cross-checks candidate sites between **PacBio and Illumina datasets**.

- Sites supported by **any Illumina sample** are excluded.
- Sites without Illumina support are retained as **high-confidence negative sites**.

---

# Pipeline for Complete Negative Control Generation

To regenerate the full negative control set:

### 1. Generate Homozygous Reference Sites
Run Sentieon Haplotyper (or GATK HaplotypeCaller) in GVCF mode:

```
1.Haplotyper_GVCF.sh
```

for both **COLO829T** and **COLO829BL**.

---

### 2. Generate Candidate Homozygous Reference BED Files

Run:

```
refhom_candidates.sh <COLO829T VCF> <COLO829BL VCF>
```

This will generate **605 BED files** under:

```
refHomCandidates/
```

---

### 3. Run the Negative Control Pipeline

The remaining steps are executed using **Snakemake**.

The pipeline requires a configuration file (example):

```
Snakemake/config/config.yaml
```

Update the config file with appropriate paths for:

- Reference genome  
  - `ref_fa`

- Path to the negative control scripts  
  - `3.negative_control`

- PacBio BAM files  
  - `pacbio[tumor_bam]`  
  - `pacbio[normal_bam]`

- Illumina BAM files  
  - `illumina[tumor_bam]`  
  - `illumina[normal_bam]`

- Pre-generated homozygous reference candidate BED files  
  - `beds[homref_candidates]`

---

### 4. Execute the Pipeline

```
snakemake --snakefile Snakemake/Snakefile --configfile Snakemake/config/config.yaml
```

---

# Output

The final **Negative Control Set** will be generated under:

```
NegativeControlSet/
```

The directory will contain BED files:

```
Step1*.bed
Step3*.bed
```

representing the filtered homozygous reference positions included in the negative control set.