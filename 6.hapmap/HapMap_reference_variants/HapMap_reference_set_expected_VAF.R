#!/usr/bin/env Rscript

#####################################################################################################################################
### HapMap Reference SNV set generation 
### with calculation of expected Variant Allele Frequencies based on mixture ratios of cell lines in SMaHT Benchmarking experiment
### by Julia Markowski
#####################################################################################################################################

# ---- Usage -------------------------------------------------------------------
# Rscript HapMap_reference_set_expected_VAF.R \
# --input_dir /path/to/vcfs \
# --output_dir /path/to/results 


# ---- Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(vcfR)
  library(dplyr)
  library(tidyr)
})

# ---- 1. Define command-line options ------------------------------------------
option_list <- list(
  make_option(
    c("-i", "--input_dir"),
    type = "character",
    help = "Input directory containing VCF files."
  ),
  make_option(
    c("-o", "--output_dir"),
    type = "character",
    help = "Output directory for results (default: same as input).",
    default = NULL
  )
)

# ---- 2. Parse arguments ------------------------------------------------------
opt <- parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options]"))

if (is.null(opt$input_dir)) {
  stop("Error: --input_dir must be provided.\nUse --help for details.")
}

input_dir <- normalizePath(opt$input_dir)

if (!dir.exists(input_dir)) {
  stop("Error: input_dir does not exist: ", input_dir)
}

output_dir <- if (!is.null(opt$output_dir)) normalizePath(opt$output_dir) else input_dir

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---- 3. Define file paths ----------------------------------------------------
HG002_vcf_path <- file.path(input_dir, "HG002_GRCh38_1_22_v4.2.1_benchmark.filtered.vcf.gz")
HG005_vcf_path <- file.path(input_dir, "HG005_GRCh38_1_22_v4.2.1_benchmark.filtered.vcf.gz")
HPRC_vcf_path <- file.path(input_dir, "hprc-v1.0-mc.grch38.vcfbub.a100k.wave.filtered.vcf.gz")

output_vcf <- file.path(output_dir, "HapMapMixture_variant_reference_set.vcf")

# ---- 4. Read VCFs  -----------------------------------------------------------
HG002_vcf <- read.vcfR(HG002_vcf_path, verbose = FALSE)
HG005_vcf <- read.vcfR(HG005_vcf_path, verbose = FALSE)
HPRC_vcf <- read.vcfR(HPRC_vcf_path, verbose = FALSE)

# ---- 5. Extract variant information and GT  ----------------------------------
HPRC_df = cbind(as.data.frame(getFIX(HPRC_vcf)),as.data.frame(extract.gt(HPRC_vcf, IDtoRowNames = FALSE)))
HG002_df = cbind(as.data.frame(getFIX(HG002_vcf)),as.data.frame(extract.gt(HG002_vcf, IDtoRowNames = FALSE)))
HG005_df = cbind(as.data.frame(getFIX(HG005_vcf)),as.data.frame(extract.gt(HG005_vcf, IDtoRowNames = FALSE)))

# GT ".|." is marked as NA, replace it with .|.
HPRC_df[, (ncol(HPRC_df) - 3):ncol(HPRC_df)] <- lapply(HPRC_df[, (ncol(HPRC_df) - 3):ncol(HPRC_df)], function(x) { x[is.na(x)] <- ".|."; x })

# ---- 6. Merge cell lines  ----------------------------------------------------
# keep all variants in merge
all_cell_lines_df <- HG002_df %>%
  full_join(HG005_df[,c("CHROM","POS","REF","ALT","HG005")], by = c("CHROM","POS","REF","ALT")) %>%
  full_join(HPRC_df[,c("CHROM","POS","REF","ALT","HG00438","HG02257","HG02486","HG02622")], by = c("CHROM","POS","REF","ALT"))

# replace missing GTs with "0|0" 
all_cell_lines_df[8:13] <- lapply(all_cell_lines_df[8:13], function(x) { x[is.na(x)] <- "0|0"; x })

# set all ID, QUAL, and FILTER to "."
all_cell_lines_df$ID = "."
all_cell_lines_df$QUAL = "."
all_cell_lines_df$FILTER = "."

# ---- 7. add expected VAF in mixture  ------------------------------------------
# add VAF per cell line, depending on homozygous REF, heterozygous, or homozygous ALT GT

assign_rawVAF <- function(gt) {
  sapply(gt, function(x) {
    # Treat NA as 0
    if (is.na(x) || x == ".|." || x == "./.") return(0)
    # Split by | or /
    alleles <- unlist(strsplit(x, split = "[|/]"))
    # Count ALT alleles (1)
    n_alt <- sum(alleles == "1")
    # Assign VAF
    if (n_alt == 2) return(1)
    if (n_alt == 1) return(0.5)
    return(0)
  })
}

all_cell_lines_df[paste0(names(all_cell_lines_df)[(ncol(all_cell_lines_df)-5):ncol(all_cell_lines_df)], "_rawVAF")] <-
  lapply(all_cell_lines_df[(ncol(all_cell_lines_df)-5):ncol(all_cell_lines_df)], assign_rawVAF)

### add adjusted VAF per cell depending on ratio in the mixture
# 0.5%: HG00438 (F, EAS)
# 2%: HG002 (M, EUR/ASH JEW)
# 2%: HG02257 (F, AFR)
# 2%: HG02486 (M, AFR)
# 10%: HG02622 (F, AFR)
# 83.5%: HG005 (M, CHB)

all_cell_lines_df$HG00438_mixVAF = all_cell_lines_df$HG00438_rawVAF*0.005
all_cell_lines_df$HG002_mixVAF = all_cell_lines_df$HG002_rawVAF*0.02
all_cell_lines_df$HG02257_mixVAF = all_cell_lines_df$HG02257_rawVAF*0.02
all_cell_lines_df$HG02486_mixVAF = all_cell_lines_df$HG02486_rawVAF*0.02
all_cell_lines_df$HG02622_mixVAF = all_cell_lines_df$HG02622_rawVAF*0.10
all_cell_lines_df$HG005_mixVAF = all_cell_lines_df$HG005_rawVAF*0.835

# add summarized VAF over all cell lines in the mixture
mixVAF_cols <- grep("_mixVAF", names(all_cell_lines_df), value = TRUE)
all_cell_lines_df$expected_mixVAF <- rowSums(all_cell_lines_df[, mixVAF_cols], na.rm = TRUE)

# ---- 8. remove variants present in background cell line HG005  ---------------
all_cell_lines_df_somatic = all_cell_lines_df[which(all_cell_lines_df$HG005_rawVAF==0),]

# ---- 9. write out VCF  -------------------------------------------------------
gt_cols <- c("HG002", "HG005", "HG00438", "HG02257", "HG02486", "HG02622")

all_cell_lines_df_somatic <- all_cell_lines_df_somatic %>%
  mutate(across(all_of(gt_cols), ~ gsub("/", "|", as.character(.))))

### Generate VAF Header
vcf_meta <- c(
  "##fileformat=VCFv4.2",
  "##INFO=<ID=HG00438,Number=1,Type=String,Description=\"Genotype in HG00438\">",
  "##INFO=<ID=HG002,Number=1,Type=String,Description=\"Genotype in HG002\">",
  "##INFO=<ID=HG02257,Number=1,Type=String,Description=\"Genotype in HG02257\">",
  "##INFO=<ID=HG02486,Number=1,Type=String,Description=\"Genotype in HG02486\">",
  "##INFO=<ID=HG02622,Number=1,Type=String,Description=\"Genotype in HG02622\">",
  "##INFO=<ID=HG005,Number=1,Type=String,Description=\"Genotype in HG005\">",
  "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Expected variant allele frequency in HapMap cell line mixture\">"
)

### Generate VAF Body
# Create INFO field with GTs and build the VCF body
vcf_body <- all_cell_lines_df_somatic %>%
  mutate(
    # build the INFO string with semicolon-separated GT info
    INFO = paste0(
      "HG00438=", HG00438, ";",
      "HG002=", HG002, ";",
      "HG02257=", HG02257, ";",
      "HG02486=", HG02486, ";",
      "HG02622=", HG02622, ";",
      "HG005=", HG005
    ),
    ID = ".",
    QUAL = ".",
    FILTER =  ".",
    FORMAT = "AF",
    HapMapMixture = expected_mixVAF
  ) %>%
  select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, HapMapMixture)

## sort variants by CHROM and POS
vcf_body$POS <- as.numeric(vcf_body$POS)
chrom_order <- sub("chr", "", vcf_body$CHROM)
chrom_order <- as.numeric(chrom_order)

vcf_body_sorted <- vcf_body[order(chrom_order, vcf_body$POS), ]

### write out VCF
options(scipen = 999) # avoid scientific notation
write_vcf <- function(vcf_meta, vcf_body, file) {
  con <- file(file, "w")
  writeLines(vcf_meta, con)
  header <- paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "HapMapMixture", sep = "\t")
  writeLines(header, con)
  write.table(vcf_body, con, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  close(con)
}

write_vcf(vcf_meta, vcf_body_sorted, output_vcf)

