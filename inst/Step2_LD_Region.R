#--------------------------------------------------------------------
# LD Matrix and Eigen-Decomposition Construction for LDER Analysis
#
# This script constructs linkage disequilibrium (LD) matrices and
# performs eigen decomposition for specified LD blocks, preparing
# inputs suitable for Linkage Disequilibrium Eigenvalue Regression (LDER).
#
# Required Input Formats:
# 1. Genotype files:
#    - PLINK binary files (.bed, .bim, .fam) organized by chromosome
#    - All files placed in a single directory (geno)
#    - Filename format example: EUR_chr_1_aligned.bed, EUR_chr_1_aligned.bim, EUR_chr_1_aligned.fam
#
# 2. GWAS aligned summary statistics files:
#    - Columns required: CHROM (chromosome), POS (position), SNP (SNP identifier),
#      ALLELE0 (reference allele), ALLELE1 (alternate allele), Z (Z-score)
#    - Organized by trait and ancestry
#    - Filename example: LDL_EUR_aligned.txt
#
# 3. LD block file:
#    - BED formatted file specifying LD blocks (chromosome, start, end positions)
#
# Input Parameter Descriptions:
# - input_dir: Main directory containing genotype and GWAS files
# - ancestry_1: Label identifying the reference ancestry (e.g., EUR)
# - ancestry_2: Label identifying the target ancestry (e.g., EAS)
# - trait: Trait being analyzed (e.g., LDL)
# - block_index: Index of the LD block to process
# - ld_block_file: File containing LD block definitions
# - ld_matrix_file: (Optional) Path to precomputed LD matrix file to skip recomputation
#
# Main functionalities:
# - Extract SNPs within specified LD blocks from genotype data
# - Compute LD matrices for both ancestries
# - Perform eigen decomposition of LD matrices
# - Generate LDER-compatible inputs from eigenvectors
# - Outputs eigen decomposition data and GWAS summary statistics for specified LD blocks
#
# Dependencies:
# - R libraries: optparse, data.table, dplyr, snpStats
#
# Example usage:
# Rscript Step4_LD_Region.R \
#   --input_dir path/to/main_directory \
#   --ancestry_1 EUR \
#   --ancestry_2 EAS \
#   --trait LDL \
#   --block_index 2 \
#   --ld_block_file path/to/grch37.eur.eas.loci.bed \
#   --ld_matrix_file NULL
#--------------------------------------------------------------------
library(optparse)
library(data.table)
library(dplyr)
library(snpStats)

# -------------------------------
# Input Parameters Specification
# -------------------------------
option_list <- list(
  make_option("--input_dir", type="character", help="Directory containing PLINK files"),

  make_option("--ancestry_1", type="character", help="Reference ancestry (e.g., EUR)"),
  make_option("--ancestry_2", type="character", help="Target ancestry (e.g., EAS)"),
  make_option("--trait", type="character", help="Trait label (e.g., LDL)"),

  make_option("--block_index", type="numeric", help="Index of the LD block to process"),
  make_option("--ld_block_file", type="character", help="File path for LD Block information"),
  make_option("--ld_matrix_file", type="character", default=NULL, help="File path for existing LD matrix (optional)")

)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Assign variables
input_dir <- opt$input_dir


ancestry_1 <- opt$ancestry_1
ancestry_2 <- opt$ancestry_2
trait <- opt$trait
block_index <- opt$block_index

dir.create(input_dir, recursive=TRUE, showWarnings=FALSE)
# Geno directory save the genotype file
geno_dir <- file.path(input_dir, "geno")
dir.create(geno_dir, recursive = TRUE, showWarnings = FALSE)
# Trait directory save the gwas file
trait_dir <- file.path(input_dir, trait)
dir.create(trait_dir , recursive = TRUE, showWarnings = FALSE)


ld_block_file <- opt$ld_block_file
ld_matrix_file <- opt$ld_matrix_file
plink_path <- opt$plink_path

# -------------------------------
# Split Region
# -------------------------------
cat("Loading LD block information...\n")
ld_block_info<-fread(ld_block_file)
chr = ld_block_info$chrom[block_index]
start = ld_block_info$start[block_index]
end = ld_block_info$stop[block_index]
cat(sprintf("Processing block index %d (chr%s:%d-%d)\n", block_index, chr, start, end))

# -------------------------------
# Read in GWAS & Reference/Match SNP
# -------------------------------
gwas_1<-fread(file.path(trait_dir, paste0(trait, "_", ancestry_1, "_aligned.txt")))
gwas_2<-fread(file.path(trait_dir, paste0(trait, "_", ancestry_2, "_aligned.txt")))

gwas_1_block<-gwas_1%>%filter(CHROM == chr)%>%filter(POS>start&POS<end)%>%arrange(as.numeric(CHROM), as.numeric(POS))
gwas_2_block<-gwas_2%>%filter(CHROM == chr)%>%filter(POS>start&POS<end)%>%arrange(as.numeric(CHROM), as.numeric(POS))

SNP_in_block = gwas_1%>%filter(CHROM == chr)%>%filter(POS>start&POS<end)%>%pull(SNP)

# ----------------------------------------------------------------------------------
# Construct LD matrix, Eigen-decomposition (Check if LD/Eigen-decomposition already finished)
# ----------------------------------------------------------------------------------
if(is.null(ld_matrix_file)){
	cat("Calculating LD matrix and eigen decomposition for ancestry 1...\n")
	ancestry_1_geno<-read.plink(file.path(geno_dir ,paste0(ancestry_1,"_chr_",chr ,"_aligned.bed")),select.snps =SNP_in_block )
	ancestry_1_geno_mat<-as(ancestry_1_geno$genotypes,"numeric")

	ancestry_1_geno_mat<-apply(ancestry_1_geno_mat,2,function(x){
		x[is.na(x)]=mean(x,na.rm=T)
		x = x - mean(x)
		x = scale(x)
		return(x)
	})
	ancestry_1_cov<-cov(ancestry_1_geno_mat)
	ancestry_1_cov_eigen <- eigen(ancestry_1_cov,symmetric = T)
	colnames(ancestry_1_cov) <- ancestry_1_geno$map$snp.name

	cat("Calculating LD matrix and eigen decomposition for ancestry 2...\n")
	ancestry_2_geno<-read.plink(file.path(geno_dir ,paste0(ancestry_2,"_chr_",chr ,"_aligned.bed")),select.snps =SNP_in_block)
	ancestry_2_geno_mat<-as(ancestry_2_geno$genotypes,"numeric")

	ancestry_2_geno_mat<-apply(ancestry_2_geno_mat,2,function(x){
		x[is.na(x)]=mean(x,na.rm=T)
		x = x - mean(x)
		x = scale(x)
		return(x)
	})
	ancestry_2_cov<-cov(ancestry_2_geno_mat)
	ancestry_2_cov_eigen <- eigen(ancestry_2_cov,symmetric = T)
	colnames(ancestry_2_cov) <- ancestry_2_geno$map$snp.name


	LD_ref_dir <- file.path(trait_dir, "LD_ref")
	if (!dir.exists(LD_ref_dir)) {
	dir.create(LD_ref_dir, recursive = TRUE, showWarnings = FALSE)
	}
	LD_eigen_out<- file.path(trait_dir,paste0("/LD_ref/LD_",block_index,".RData"))
	save(ancestry_1_cov,ancestry_1_cov_eigen,ancestry_2_cov,ancestry_2_cov_eigen ,file = LD_eigen_out)
	cat("LD matrix and eigen decomposition saved.\n")

}else{
	cat("Loading existing LD matrix and eigen decomposition...\n")
	loaded_vars<-load(ld_matrix_file)
	required_vars <- c("ancestry_1_cov", "ancestry_1_cov_eigen", "ancestry_2_cov", "ancestry_2_cov_eigen")
  	missing_vars <- setdiff(required_vars, loaded_vars)
	if (length(missing_vars) > 0) {
		stop(paste("Missing required variables in LD matrix file:", paste(missing_vars, collapse=", ")))
	}

}

# ----------------------------------------------------
# 				Write out GWAS files in blocks
# ----------------------------------------------------

# Write out gwas block
plink_snp_order <- colnames(ancestry_1_cov)
gwas_1_block <- gwas_1_block %>%
  arrange(match(SNP, plink_snp_order))
gwas_2_block <- gwas_2_block %>%
  arrange(match(SNP, plink_snp_order))

stopifnot(identical(gwas_1_block$SNP, plink_snp_order))
stopifnot(identical(gwas_2_block$SNP, plink_snp_order))

GWAS_block_dir <- file.path(trait_dir, paste0(trait,"_Block"))
if (!dir.exists(GWAS_block_dir)) {
  dir.create(GWAS_block_dir, recursive = TRUE, showWarnings = FALSE)
}

fwrite(gwas_1_block,file.path(GWAS_block_dir, paste0(trait, "_", ancestry_1, "_block_",block_index,".txt")))
fwrite(gwas_2_block,file.path(GWAS_block_dir, paste0(trait, "_", ancestry_2, "_block_",block_index,".txt")))


# ----------------------------------------------------
# Construct Input for LDER by eigenvectors of LD matrix
# ----------------------------------------------------
lder_preprocess<-function(temp,z,eff.num){
        temp$values[temp$values<1e-6] <- 0
        U <- temp$vectors
        V <- diag(temp$values)
        V.inv <-  1/V
        V.inv[which(V.inv==Inf)] <- 0
        if(is.null(eff.num)){eff.num <- 1}
        eigen.mat<- sqrt(V.inv[1:eff.num,1:eff.num])%*%(t(U)[1:eff.num,])
        lam <- diag(V)[1:eff.num] ## export
        z <- z
        x <- eigen.mat%*%z ## export
		return(list(x=x,z=z,lam=lam))
}

num_SNP = nrow(gwas_1_block)
LDER_ancestry_1<-lder_preprocess(ancestry_1_cov_eigen,gwas_1_block$Z,num_SNP)
LDER_ancestry_2<-lder_preprocess(ancestry_2_cov_eigen,gwas_2_block$Z,num_SNP)

LDER_out<-file.path(GWAS_block_dir ,paste0(trait,"_LDER_block_",block_index,".RData"))
save(LDER_ancestry_1,LDER_ancestry_2,file = LDER_out)





#ancestry_1_geno<-read.plink(file.path(geno_dir ,paste0(ancestry_1,"_chr_",chr ,"_aligned.bed")),select.snps =SNP_in_block )
#ancestry_2_geno<-read.plink(file.path(geno_dir ,paste0(ancestry_2,"_chr_",chr ,"_aligned.bed")),select.snps =SNP_in_block)
# sum(ancestry_1_geno$map$allele.1==gwas_1_block$ALLELE0)
# sum(ancestry_2_geno$map$allele.1==gwas_2_block$ALLELE0)

#input_dir<-"/net/fantasia/home/borang/MALGC/pipeline_example"
#ancestry_1<- "EUR"
#ancestry_2<- "EAS"
#trait<-"LDL"
#output_dir<-"/net/fantasia/home/borang/MALGC/pipeline_example"
#ld_block_file<-"/net/fantasia/home/borang/MALGC/ld_blocks/grch37.eur.eas.loci.bed"
# block_index<-1
# Rscript ~/MALGC/MALGC_software/Data_Process/Step4_LD_Region.R --input_dir /net/fantasia/home/borang/MALGC/pipeline_example --output_dir /net/fantasia/home/borang/MALGC/pipeline_example --ancestry_1 EUR --ancestry_2 EAS --trait LDL --block_index 2 --ld_block_file /net/fantasia/home/borang/MALGC/ld_blocks/grch37.eur.eas.loci.bed 