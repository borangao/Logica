#--------------------------------------------------------------------
# LDER Intercept Estimation Script
#
# This script performs Linkage Disequilibrium Eigenvalue Regression
# (LDER) analysis using GWAS summary statistics for two specified ancestries.
# It aggregates LD block data, computes intercepts using the 'lder' function
# from the LDER R package, and saves intercept estimates for downstream
# genetic correlation analyses.
#
# Required Input Formats:
# - Preprocessed LDER GWAS blocks:
#   - RData (.RData) files generated from prior eigen decomposition analyses.
#   - Files stored in the directory structure: input_dir/[trait]/[trait]_Block/
#   - Filename example: LDL_LDER_block_1.RData
#
# Input Parameter Descriptions:
# - input_dir: Directory containing trait-specific subdirectories with GWAS blocks
# - ancestry_1: Label identifying the reference ancestry (e.g., EUR)
# - n_ancestry_1: Sample size for the first ancestry
# - ancestry_2: Label identifying the target ancestry (e.g., EAS)
# - n_ancestry_2: Sample size for the second ancestry
# - trait: Trait being analyzed (e.g., LDL)
#
# Main functionalities:
# 1. Load LD block data (.RData files) for each ancestry.
# 2. Perform LDER analysis to estimate intercept values.
# 3. Save intercept estimates to an RData file for downstream use.
#
# Dependencies:
# - R libraries: optparse, data.table, dplyr, LDER
# - Installation from GitHub: devtools::install_github('borangao/LDER')
#
# Example usage:
# Rscript Step3_LDER.R \
#   --input_dir /path/to/input \
#   --ancestry_1 EUR --n_ancestry_1 343621 \
#   --ancestry_2 EAS --n_ancestry_2 237613 \
#   --trait LDL
#
# Output:
# - Intercept estimates saved as '[Trait]_intercept.RData' in input_dir/[Trait]/
# ----------------------------------------------------------------------------
library(optparse)
library(data.table)
library(dplyr)

devtools::install_github('borangao/LDER')
library(LDER)
# -------------------------------
# Input Parameters Specification
# -------------------------------
option_list <- list(
  make_option("--input_dir", type="character", help="Directory containing PLINK files"),
  make_option("--ancestry_1", type="character", help="Reference ancestry (e.g., EUR)"),
  make_option("--n_ancestry_1", type="numeric", help="N Sample of first ancestry"),
  make_option("--ancestry_2", type="character", help="Target ancestry (e.g., EAS)"),
  make_option("--n_ancestry_2", type="numeric", help="N Sample of second ancestry"),
  make_option("--trait", type="character", help="Trait label (e.g., LDL)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

input_dir <- opt$input_dir

ancestry_1 <- opt$ancestry_1
ancestry_2 <- opt$ancestry_2

n_ancestry_1 <- opt$n_ancestry_1
n_ancestry_2 <- opt$n_ancestry_2

trait <- opt$trait

GWAS_block_dir<-file.path(input_dir, trait, paste0(trait, "_Block"))

setwd(GWAS_block_dir)
cat("Loading LD block data from directory:", GWAS_block_dir, "\n")

LDER_list <- list.files(path = GWAS_block_dir, 
                        pattern = paste0(trait, "_LDER_block_.*\\.RData$"),
                        full.names = TRUE)

LDER_ancestry_1_all<-list()
LDER_ancestry_2_all<-list()

for(LD_BLOCK_NUM in 1:length(LDER_list)){
  LDER_out<-LDER_list[LD_BLOCK_NUM]
  load(LDER_out)
  
  LDER_ancestry_1_all[[LD_BLOCK_NUM ]]<-LDER_ancestry_1
  LDER_ancestry_2_all[[LD_BLOCK_NUM ]]<-LDER_ancestry_2
}

set.seed(1103)
res_lder_ancestry_1 <- lder(stats=LDER_ancestry_1_all,n.gwas=n_ancestry_1 ,cores=20,a=NULL,rough=T,twostage=FALSE,type='none')
cat("Completed analysis for", ancestry_1, "Intercept estimate:", res_lder_ancestry_1$inf, "\n")

res_lder_ancestry_2 <- lder(stats=LDER_ancestry_2_all,n.gwas=n_ancestry_2 ,cores=20,a=NULL,rough=T,twostage=FALSE,type='none')
cat("Completed analysis for", ancestry_2, "Intercept estimate:", res_lder_ancestry_2$inf, "\n")


ancestry_1_lder_intercept<-res_lder_ancestry_1$inf
ancestry_2_lder_intercept<-res_lder_ancestry_2$inf

save(
ancestry_1_lder_intercept,
ancestry_2_lder_intercept,
file= file.path(input_dir, trait,paste0(trait, "_intercept.RData")))
cat("Intercept estimates successfully saved to", file.path(GWAS_block_dir,paste0(trait, "_intercept.RData")), "\n")


# Rscript /net/fantasia/home/borang/MALGC/MALGC_software/Data_Process/Step3_LDER.R --input_dir /net/fantasia/home/borang/MALGC/pipeline_example  --ancestry_1 EUR --n_ancestry_1 343621 --ancestry_2 EAS --n_ancestry_2 237613  --trait LDL
