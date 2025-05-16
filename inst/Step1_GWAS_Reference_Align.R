#--------------------------------------------------------------------
# GWAS and Reference Panel QC and SNP Alignment Script
#
# This script performs quality control (QC) on GWAS summary statistics
# and genotype reference panel data for two specified ancestries. 
# It aligns SNPs between the GWAS datasets and corresponding reference
# panels, filters for common SNPs, and prepares genotype files suitable
# for downstream analyses.
#
# Required Input Formats:
# 1. Reference genotype files:
#    - PLINK binary files (.bed, .bim, .fam) organized by chromosome
#    - All files placed in a single directory per ancestry
#    - Filename format example: EUR_chr_1.bed, EUR_chr_1.bim, EUR_chr_1.fam
#
# 2. GWAS summary statistics files:
#    - Required columns: CHROM (chromosome), POS (position), SNP (SNP identifier),
#      ALLELE0 (reference allele), ALLELE1 (alternate allele), Z (Z-score), N (sample size)
#    - Optional columns: BETA (effect size), SE (standard error), A1FREQUENCY (allele frequency)
#
# Input Parameter Descriptions:
# - input_dir_1: Directory path containing PLINK files for ancestry 1.
# - input_prefix_1: Prefix for PLINK file names for ancestry 1 (e.g., EUR_chr).
# - ancestry_1: Label identifying the first ancestry (e.g., EUR).
# - gwas_1: Path to GWAS summary statistics file for the first ancestry.
# - input_dir_2: Directory path containing PLINK files for ancestry 2.
# - input_prefix_2: Prefix for PLINK file names for ancestry 2 (e.g., EAS_chr).
# - ancestry_2: Label identifying the second ancestry (e.g., EAS).
# - gwas_2: Path to GWAS summary statistics file for the second ancestry.
# - trait: Name of the trait being analyzed (e.g., LDL).
# - output_dir: Directory where the output files will be stored.
# - plink_path: Full path to the PLINK 2 executable.
# - skip_geno_qc: Boolean flag indicating whether to skip genotype QC steps (default FALSE).
#
# Main functionalities:
# - QC on GWAS summary statistics (filters ambiguous SNPs, missing data,
#   allele frequency, and MHC region)
# - QC on reference genotype data using PLINK (filters MAF, genotype rate,
#   HWE)
# - Alignment of SNPs across two ancestries and GWAS summary statistics
# - Generates standardized outputs:
#    - QC'd genotype files (PLINK binary format by chromosome)
#    - Aligned GWAS summary statistics
#    - SNP lists for further analysis
#
# Dependencies:
# - R libraries: optparse, data.table, dplyr
# - External software: PLINK 2
#
# Example usage:
# Rscript Step1_GWAS_Reference_Align.R \
#   --input_dir_1 path/to/ancestry1 \
#   --input_prefix_1 EUR_chr \
#   --ancestry_1 EUR \
#   --gwas_1 EUR_GWAS.txt \
#   --input_dir_2 path/to/ancestry2 \
#   --input_prefix_2 EAS_chr \
#   --ancestry_2 EAS \
#   --gwas_2 EAS_GWAS.txt \
#   --trait LDL \
#   --output_dir path/to/output \
#   --plink_path /path/to/plink2 \
#   --skip_geno_qc FALSE
# Rscript Step1_GWAS_Reference_Align.R --input_dir_1 /net/fantasia/home/borang/MALGC/ukb_bbj_ref/EUR/ --input_prefix_1 EUR_chr --ancestry_1 EUR  --gwas_1 /net/fantasia/home/borang/MALGC/real_data/European/UKBB/UKBB_QC/LDL_common.txt --input_dir_2 /net/fantasia/home/borang/MALGC/ukb_bbj_ref/EAS/ --input_prefix_2 EAS_chr --ancestry_2 EAS --gwas_2 /net/fantasia/home/borang/MALGC/real_data/Asian/Meta/LDL_common.txt --trait LDL  --output_dir /net/fantasia/home/borang/MALGC/pipeline_example --plink_path /net/fantasia/home/borang/software/plink2 --skip_geno_qc TRUE
#--------------------------------------------------------------------

library(optparse)
library(data.table)
library(dplyr)
## This script performs quality control (QC) on reference panel
## Required plink file for both ancestries save in by chromosome format
## Output plink file for each ancestry in [output_dir]/geno/ director with plink file saved in by chromosome format

# 1. Specificy Input
option_list <- list(
  make_option("--input_dir_1", type="character", help="Directory with PLINK files of the first ancestry"),
  make_option("--input_prefix_1", type="character", help="Prefix for PLINK files (e.g., EUR_chr)"),
  make_option("--ancestry_1", type="character", help="Ancestry label 1 (e.g., EUR)"),
  make_option("--gwas_1", type="character", help="Input GWAS summary stats file", metavar="character"),

  make_option("--input_dir_2", type="character", help="Directory with PLINK files of the second ancestry"),
  make_option("--input_prefix_2", type="character", help="Prefix for PLINK files (e.g., EAS_chr)"),
  make_option("--ancestry_2", type="character", help="Ancestry label 2 (e.g., EAS)"),
  make_option("--gwas_2", type="character", help="Input GWAS summary stats file", metavar="character"),
   
  make_option("--trait", type="character", help="Trait label", metavar="character"),

  make_option("--output_dir", type="character", help="Output directory"),
  make_option("--plink_path", type="character", help="Full path to user-specified PLINK binary"),
  make_option("--skip_geno_qc", action="store_true", default=FALSE,
              help="Set this flag to skip genotype QC if genotype files already exist.")

)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

input_dir_1 <- opt$input_dir_1
input_prefix_1 <- opt$input_prefix_1
ancestry_1 <- opt$ancestry_1

input_dir_2 <- opt$input_dir_2
input_prefix_2 <- opt$input_prefix_2
ancestry_2 <- opt$ancestry_2

trait <- opt$trait

output_dir <- opt$output_dir
plink_path<-opt$plink_path

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
# Geno directory save the genotype file
geno_dir <- file.path(output_dir, "geno")
dir.create(geno_dir, recursive = TRUE, showWarnings = FALSE)
# Trait directory save the gwas file
trait_dir <- file.path(output_dir, trait)
dir.create(trait_dir , recursive = TRUE, showWarnings = FALSE)



# QC thresholds
maf_thresh <- 0.05
geno_thresh <- 0.05
hwe_thresh <- 1e-7

# 2. List of SNPs to Keep
generate_snp_list <- function(input_dir, input_prefix) {
  bim_all <- lapply(1:22, function(chr) {
    bim_file <- file.path(input_dir, paste0(input_prefix, "_", chr, ".bim"))
    fread(bim_file)
  }) %>% rbindlist()

  setnames(bim_all, c("CHROM", "SNP", "CM", "POS", "ALLELE1", "ALLELE0"))

  ambiguous_pairs <- c("AT", "TA", "CG", "GC")
  bim_all[, allele_pair := paste0(ALLELE1, ALLELE0)]
  ambiguous_snps <- bim_all[allele_pair %in% ambiguous_pairs, SNP]

  bim_filtered <- bim_all[!SNP %in% ambiguous_snps] %>%
    .[!(CHROM == 6 & POS >= 25000000 & POS <= 36000000)] %>%
    .[!duplicated(paste(CHROM, POS))]

  return(unique(bim_filtered$SNP))
}
if (!opt$skip_geno_qc) {

cat("Generating SNP keep list...\n")
snp_keep_1 <- generate_snp_list(input_dir_1, input_prefix_1)
snp_keep_2 <- generate_snp_list(input_dir_2, input_prefix_2)

snp_keep_file_1 <- file.path(geno_dir, paste0(ancestry_1, "_snps_to_keep.txt"))
snp_keep_file_2 <- file.path(geno_dir, paste0(ancestry_2, "_snps_to_keep.txt"))

fwrite(data.table(SNP=snp_keep_1), snp_keep_file_1, col.names=FALSE)
fwrite(data.table(SNP=snp_keep_2), snp_keep_file_2, col.names=FALSE)

# 3. QC on PLINK Files

ancestries <- list(
  list(input_dir = input_dir_1, input_prefix = input_prefix_1, ancestry = ancestry_1, snp_keep_file = snp_keep_file_1),
  list(input_dir = input_dir_2, input_prefix = input_prefix_2, ancestry = ancestry_2, snp_keep_file = snp_keep_file_2)
)

for (anc in ancestries) {
  input_dir <- anc$input_dir
  input_prefix <- anc$input_prefix
  ancestry <- anc$ancestry
  snp_keep_file <- anc$snp_keep_file

  for (chr in 1:22) {
    input_plink_prefix <- file.path(input_dir, paste0(input_prefix, "_", chr))
    output_plink_prefix <- file.path(geno_dir, paste0(ancestry, "_chr_", chr, "_qc"))

    plink_cmd <- sprintf(
      "%s --bfile %s --extract %s --maf %.2f --geno %.2f --hwe %.1e --make-bed --out %s",
      plink_path, input_plink_prefix, snp_keep_file, maf_thresh, geno_thresh, hwe_thresh, output_plink_prefix
    )

    cat(sprintf("\nRunning PLINK QC for ancestry %s chromosome %d:\n%s\n", ancestry, chr, plink_cmd))
    system(plink_cmd)

    # Verify output files
    bim_qc_file <- paste0(output_plink_prefix, ".bim")
    if (file.exists(bim_qc_file)) {
      bim <- fread(bim_qc_file)
      setnames(bim, c("CHROM", "SNP", "CM", "POS", "ALLELE1", "ALLELE0"))
      bim[, SNP := paste0(CHROM, "_", POS)]  # Standardize SNP names
      fwrite(bim, bim_qc_file, sep=" ", col.names=FALSE)
      cat(sprintf("Ancestry %s chromosome %d QC and SNP renaming completed.\n", ancestry, chr))
    } else {
      warning(sprintf("PLINK QC failed for ancestry %s chromosome %d. Please check manually.", ancestry, chr))
    }
  }
}
}else{
 cat("Skipping genotype QC as requested (--skip_geno_qc is TRUE)...\n")
}


## This script performs quality control (QC) on GWAS summary statistics data.
## Required columns: Chromosome (CHROM), Base-pair position (POS), SNP identifier (SNP), Reference allele (ALLELE0), Alternative allele (ALLELE1), Z-score (Z), Sample size (N).
## Optional columns: Effect size (BETA), Standard error (SE), Minor allele frequency (MAF).

## QC Steps include:
## 1. Standardize column names to uppercase.
## 2. Convert alleles (ALLELE0 and ALLELE1) to uppercase.
## 3. If SNP identifier (SNP) is missing, generate SNP IDs using CHROM and POS.
## 4. Remove SNPs with missing essential information (CHROM, POS, SNP, ALLELE0, ALLELE1).
## 5. Keep only standard nucleotide alleles (A, T, C, G).
## 6. Remove strand-ambiguous SNPs (allele pairs: A/T, T/A, C/G, G/C).
## 7. Exclude multi-allelic SNPs.
## 8. Exclude SNPs located within the MHC region on chromosome 6 (25Mb–36Mb).
## 9. If A1FREQUENCY column is present, retain only SNPs with A1FREQUENCY between 0.05 and 0.95.

perform_qc_gwas <- function(data) {
  
  required_cols <- c("CHROM", "POS", "ALLELE0", "ALLELE1", "N")
  
  colnames(data) <- toupper(colnames(data))
  
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  has_z <- "Z" %in% colnames(data)
  has_beta_se <- all(c("BETA", "SE") %in% colnames(data))
  
  if (!(has_z || has_beta_se)) {
    stop("GWAS data must contain either 'Z' or both 'BETA' and 'SE' columns.")
  }
  
  data <- data %>% mutate(SNP = paste(CHROM, POS, sep = "_"))
  
  data <- data %>%
    mutate(
      ALLELE0 = toupper(ALLELE0),
      ALLELE1 = toupper(ALLELE1)
    )
  
  strand_ambiguous_pairs <- c("AT", "TA", "CG", "GC")
  
  initial_count <- nrow(data)
  
  data_qc <- data %>%
    filter(
      !is.na(CHROM), !is.na(POS), !is.na(SNP),
      !is.na(ALLELE0), !is.na(ALLELE1),
      ALLELE0 %in% c("A", "T", "C", "G"),
      ALLELE1 %in% c("A", "T", "C", "G"),
      !(paste0(ALLELE0, ALLELE1) %in% strand_ambiguous_pairs)
    ) %>%
    add_count(CHROM, POS) %>%
    filter(n == 1) %>%
    select(-n) %>%
    filter(!(CHROM == "6" & POS >= 25000000 & POS <= 36000000))
  
    # USE A1FREQUENCY to filter out SNPs
  if ("A1FREQUENCY" %in% colnames(data_qc)) {
    data_qc <- data_qc %>% filter(A1FREQUENCY >= 0.05 & A1FREQUENCY <= 0.95)
  }
  
  if (!has_z && has_beta_se) {
    data_qc <- data_qc %>% mutate(Z = BETA / SE)
    cat("Z column computed from BETA and SE.\n")
  }
  
  final_count <- nrow(data_qc)
  cat(sprintf("Retained %d SNPs out of %d after QC.\n", final_count, initial_count))
  
  return(data_qc)
}


# read in GWAS data
gwas_1 <- fread(opt$gwas_1)
gwas_2 <- fread(opt$gwas_2)
#gwas_1<-fread("/net/fantasia/home/borang/MALGC/real_data/European/UKBB/UKBB_QC/LDL_common.txt")
#gwas_2<-fread("/net/fantasia/home/borang/MALGC/real_data/Asian/Meta/LDL_common.txt")
# perform QC
cat(paste0("Perform GWAS QC for ancestry ",ancestry_1))
gwas_1 <- perform_qc_gwas(gwas_1)
gwas_1 <-gwas_1%>%mutate(CHROM = as.numeric(CHROM))
cat(paste0("Perform GWAS QC for ancestry ",ancestry_2))
gwas_2 <- perform_qc_gwas(gwas_2)
gwas_2 <-gwas_2%>%mutate(CHROM = as.numeric(CHROM))
# write out
#fwrite(gwas_qc, paste0(opt$output,"/",opt$trait,"_",opt$ancestry,"_QC.txt"))

#######################################################################################
#
#               Match GWAS with Reference Panel SNPs
#
#
######################################################################################
# -------------------------------
# Step 1: Load and Prepare Data
# -------------------------------
geno_dir<-"/net/fantasia/home/borang/MALGC/pipeline_example/geno"
ancestry_1<-"EUR"
ancestry_2<-"EAS"
# Load and combine BIM files from ancestry 1 (reference)

cat("Loading and combining BIM files from ancestry 1 (reference)...\n")
bim_ancestry_1 <- lapply(1:22, function(chr) {
  fread(file.path(geno_dir, paste0(ancestry_1, "_chr_", chr, "_qc.bim")))
}) %>% rbindlist()
setnames(bim_ancestry_1, c("CHROM", "SNP", "CM", "POS", "ALLELE1", "ALLELE0"))

##Note here for bim file 5th column ALT ("A1" in PLINK 1.x) allele code, 6th column REF ("A2" in PLINK 1.x) allele code


# Load and combine BIM files from ancestry 2
cat("Loading and combining BIM files from ancestry 2...\n")
bim_ancestry_2 <- lapply(1:22, function(chr) {
  fread(file.path(geno_dir, paste0(ancestry_2, "_chr_", chr, "_qc.bim")))
}) %>% rbindlist()
setnames(bim_ancestry_2, c("CHROM", "SNP", "CM", "POS", "ALLELE1", "ALLELE0"))

# Load GWAS summary statistics
cat("Loading GWAS summary statistics for ancestry 1...\n")
gwas_1_info <- gwas_1 %>% select(CHROM, SNP, POS, ALLELE1, ALLELE0)
gwas_2_info <- gwas_2 %>% select(CHROM, SNP, POS, ALLELE1, ALLELE0)

# -------------------------------
# Step 2: Align SNPs across ancestries and GWAS
# -------------------------------
cat("Aligning SNPs across ancestries and GWAS...\n")
align_to_reference <- function(target_data, reference_data, target_suffix) {
  merged <- inner_join(target_data, reference_data, by = c("CHROM", "POS"),
                       suffix = c(paste0("_", target_suffix), "_ref"))
  if(nrow(merged) == 0) stop("No overlapping SNPs found.")

  flip <- ifelse(
    merged[[paste0("ALLELE0_", target_suffix)]] == merged$ALLELE0_ref &
    merged[[paste0("ALLELE1_", target_suffix)]] == merged$ALLELE1_ref, FALSE,
    ifelse(merged[[paste0("ALLELE0_", target_suffix)]] == merged$ALLELE1_ref &
           merged[[paste0("ALLELE1_", target_suffix)]] == merged$ALLELE0_ref, TRUE, NA)
  )

  aligned <- merged %>% filter(!is.na(flip))
  if(nrow(aligned) == 0) stop("No SNPs aligned after filtering.")

  reference_aligned <- aligned %>%
    select(CHROM, POS, ends_with("_ref")) %>%
    rename_with(~ gsub("_ref$", "", .x))

  return(reference_aligned)
}

# Align reference panels and GWAS data
common_snps <- align_to_reference(bim_ancestry_1, bim_ancestry_2, ancestry_1)
common_snps <- align_to_reference(common_snps, gwas_1_info, ancestry_1)
common_snps <- align_to_reference(common_snps, gwas_2_info, ancestry_1)

cat("Number of common SNPs identified:", nrow(common_snps), "\n")

# -------------------------------
# Step 3: Filter SNPs to common set
# -------------------------------
cat("Filtering SNPs to common set...\n")

snp_list <- common_snps$SNP
bim_ancestry_1_sub <- bim_ancestry_1 %>% filter(SNP %in% snp_list)
bim_ancestry_2_sub <- bim_ancestry_2 %>% filter(SNP %in% snp_list)
gwas_1_sub <- gwas_1 %>% filter(SNP %in% snp_list)
gwas_2_sub <- gwas_2 %>% filter(SNP %in% snp_list)

# -------------------------------
# Step 4: Align GWAS effect alleles
# -------------------------------
cat("Aligning GWAS effect alleles...\n")
align_to_reference_bim_GWAS <- function(gwas_data, bim_data) {
  
  merged <- inner_join(gwas_data, bim_data, by = c("CHROM", "POS"))

  if (nrow(merged) == 0) stop("No overlapping SNPs (CHROM, POS).")

  # 1. flip or not
flip <- ifelse(
  merged$ALLELE0.x == merged$ALLELE0.y & merged$ALLELE1.x == merged$ALLELE1.y, FALSE,
  ifelse(merged$ALLELE0.x == merged$ALLELE1.y & merged$ALLELE1.x == merged$ALLELE0.y, TRUE, NA)
)

merged$flip <- flip

# 2. remove ummatched snps
aligned <- merged %>% filter(!is.na(flip))

if (nrow(aligned) == 0) stop("No SNPs aligned after allele matching.")

# 3. flip alleles, Z
aligned <- aligned %>%
  mutate(
    ALLELE0 = ifelse(flip, ALLELE1.x, ALLELE0.x),
    ALLELE1 = ifelse(flip, ALLELE0.x, ALLELE1.x),
    Z = ifelse(flip, -Z, Z) 
  )

# flip beta and A1FREQ if exists
if ("BETA" %in% names(aligned)) {
  aligned <- aligned %>% mutate(BETA = ifelse(flip, -BETA, BETA))
}

if ("A1FREQ" %in% names(aligned)) {
  aligned <- aligned %>% mutate(A1FREQ = ifelse(flip, 1 - A1FREQ, A1FREQ))
}

# 4. Keep columns 
required_cols <- c("CHROM", "POS", "SNP.x", "ALLELE0", "ALLELE1", "Z","N")
optional_cols <- intersect(names(aligned), c("BETA", "SE", "pval",  "A1FREQ"))

gwas_aligned <- aligned %>%
  select(all_of(c(required_cols, optional_cols))) %>%
  rename_with(~ gsub("\\.x$", "", .x), everything())%>%
    arrange(as.numeric(CHROM), as.numeric(POS))
  return(gwas_aligned)
}



gwas_1_aligned <- align_to_reference_bim_GWAS(gwas_1_sub, bim_ancestry_1_sub)
gwas_2_aligned <- align_to_reference_bim_GWAS(gwas_2_sub, bim_ancestry_1_sub)

# -------------------------------
# Step 5: Output files
# -------------------------------
cat("Writing aligned GWAS data to output files...\n")
fwrite(gwas_1_aligned, file.path(trait_dir, paste0(trait, "_", ancestry_1, "_aligned.txt")), sep=" ")
fwrite(gwas_2_aligned, file.path(trait_dir, paste0(trait, "_", ancestry_2, "_aligned.txt")), sep=" ")

# Generate SNP reference and list files
allele_ref <- bim_ancestry_1_sub %>% select(SNP, ALLELE1, ALLELE0) %>% rename(A1 = ALLELE1, A2 = ALLELE0) ## Reference allele is the A2, which is the 3rd column
common_snp_file <- file.path(geno_dir, "common_snps.txt")
allele_ref_file <- file.path(geno_dir, "common_snps_allele_ref.txt")

fwrite(data.frame(SNP=allele_ref$SNP), common_snp_file, col.names=FALSE)
fwrite(allele_ref, allele_ref_file, col.names=FALSE,sep=" ")

# -------------------------------
# Step 6: PLINK Subset
# -------------------------------
for(chr in 1:22){
  system(sprintf("%s --bfile %s --extract %s --keep-allele-order --make-bed --out %s",
                 plink_path,
                 file.path(geno_dir, paste0(ancestry_1, "_chr_", chr, "_qc")),
                 common_snp_file,
                 file.path(geno_dir, paste0(ancestry_1, "_chr_", chr, "_aligned"))))

  system(sprintf("%s --bfile %s --extract %s --keep-allele-order --ref-allele %s 3 1  --make-bed --out %s",
                 plink_path,
                 file.path(geno_dir, paste0(ancestry_2, "_chr_", chr, "_qc")),
                 common_snp_file,
                 allele_ref_file,
                 file.path(geno_dir, paste0(ancestry_2, "_chr_", chr, "_aligned"))))
}
cat("All chromosomes have been processed successfully!\n")
