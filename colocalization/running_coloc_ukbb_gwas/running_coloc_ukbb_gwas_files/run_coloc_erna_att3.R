library("data.table")
library("coloc")
library("optparse")

option_list = list(
  make_option(c("-d", "--in_dir"), type = "character", default = NULL, help = "directory of the input file", metavar = "character"),
  make_option(c("-t", "--tissue"), type = "character", default = NULL, help = "tissue with which to run coloc", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

in_dir <- opt$in_dir
tissue <- opt$tissue

out_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/gwas_snp_coloc/att3"

pheno_df <- fread("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/gwas_snp_coloc/att2/phenotype_manifest.tsv")
pheno_df <- as.data.frame(pheno_df)
pheno_df$pheno_full <- ifelse(pheno_df$modifier != "", paste(pheno_df$phenocode, pheno_df$modifier, sep = "_"), pheno_df$phenocode)

erna_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/gwas_snp_coloc/erna_all_tissues_eqtls_matrix_eqtl.with_coord.arm.filt_fdr.txt"
erna_file <- fread(erna_path, header = TRUE, sep = "\t", quote = "")
erna_df <- as.data.frame(erna_file)
names(erna_df)[1] <- "varid"
erna_df <- erna_df[(erna_df$tissue == tissue),]

maf_erna_path <- "/home/bettimj/gamazon_rotation/gtex8_geno/liftover_hg19/plink2.afreq"
maf_erna_file <- fread(maf_erna_path, header = TRUE, sep = "\t", quote = "")
maf_erna_df <- as.data.frame(maf_erna_file)

erna_df <- merge(erna_df, maf_erna_df, by.x = "varid", by.y = "ID")

# Initialize an empty list to store results
results_list <- c()
pheno_locus_list <- c()
tissue_list <- c()

for (gwas_file_name in list.files(in_dir, pattern = "*.txt")) {
  print(gwas_file_name)
  file_name <- gwas_file_name
  file_name <- strsplit(file_name, "-")
  pheno <- unlist(lapply(file_name, `[`, 1))
  pheno_type <- ""
  
  pheno_df_filt <- pheno_df[(pheno_df$phenocode %in% pheno),]
  pheno_df_filt <- pheno_df_filt[!is.na(pheno_df_filt$phenocode),]
  pheno_df_filt <- pheno_df_filt[1,]
  if (any(is.na(pheno_df_filt$n_controls_AFR)) & 
      any(is.na(pheno_df_filt$n_controls_AMR)) & 
      any(is.na(pheno_df_filt$n_controls_CSA)) & 
      any(is.na(pheno_df_filt$n_controls_EAS)) & 
      any(is.na(pheno_df_filt$n_controls_EUR)) & 
      any(is.na(pheno_df_filt$n_controls_MID))) {
    pheno_type <- "continuous"
  } else {
    pheno_type <- "binary"
  }
  
  gwas_file <- fread(file.path(in_dir, gwas_file_name), header = TRUE, sep = "\t", quote = "")
  gwas_df <- as.data.frame(gwas_file)
  
  # Remove duplicate SNPs
  gwas_df <- gwas_df[!duplicated(gwas_df$varid),]
  
  # Convert MAF to numeric and remove NA values
  gwas_df$minor_AF <- as.numeric(gwas_df$minor_AF)
  gwas_df <- gwas_df[!is.na(gwas_df$minor_AF),]
  
  # Fix MAF issue by explicitly filtering for MAF between 0 and 1
  gwas_df <- gwas_df[(gwas_df$minor_AF > 0 & gwas_df$minor_AF < 1),]
  gwas_df <- na.omit(gwas_df)
  
  locus <- unique(gwas_df$arm_vector)
  
  eqtl_locus <- erna_df[(erna_df$arm_vector == locus),]
  eqtl_locus_tissue <- eqtl_locus[(eqtl_locus$tissue == tissue),]
  eqtl_locus_tissue <- eqtl_locus_tissue[order(eqtl_locus_tissue$pvalue),]
  eqtl_locus_tissue <- eqtl_locus_tissue[!duplicated(eqtl_locus_tissue$varid),]
  merged_df <- merge(gwas_df, eqtl_locus_tissue, by = "varid")
  merged_df <- merged_df[!duplicated(merged_df$varid),]
  
  if (nrow(merged_df) > 0) {
    if (pheno_type == "continuous") {
      print(pheno_type)
      dataset1 <- list(snp = merged_df$varid, beta = merged_df$beta.x, varbeta = merged_df$se^2, pvalues = merged_df$pval, type = "quant", N = nrow(gwas_df), MAF = merged_df$minor_AF)
      dataset2 <- list(snp = merged_df$varid, beta = merged_df$beta.y, pvalues = merged_df$pvalue, type = "quant", N = nrow(eqtl_locus_tissue), MAF = merged_df$ALT_FREQS)
      print(length(merged_df$varid))
    } else if (pheno_type == "binary") {
      print(pheno_type)
      n_cases <- pheno_df_filt$n_cases_full_cohort_both_sexes
      n_controls <- rowSums(pheno_df_filt[, c("n_controls_AFR", "n_controls_AMR", "n_controls_CSA", "n_controls_EAS", "n_controls_EUR", "n_controls_MID")], na.rm = TRUE)
      print(n_cases)
      print(n_controls)
      prop_cases <- n_cases / (n_cases + n_controls)
      print(prop_cases)
      dataset1 <- list(snp = merged_df$varid, beta = merged_df$beta.x, varbeta = merged_df$se^2, pvalues = merged_df$pval, type = "cc", s = prop_cases, N = nrow(gwas_df), MAF = merged_df$minor_AF)
      dataset2 <- list(snp = merged_df$varid, beta = merged_df$beta.y, pvalues = merged_df$pvalue, type = "quant", N = nrow(eqtl_locus_tissue), MAF = merged_df$ALT_FREQS)
      print(length(merged_df$varid))
    }
    tryCatch({
      result <- coloc.abf(dataset1, dataset2)
      Posterior_Probability = result$summary[6]
      results_list <- c(results_list, Posterior_Probability)
      pheno_locus <- paste(pheno, locus, sep = "_")
      pheno_locus_list <- c(pheno_locus_list, pheno_locus)
      tissue_list <- c(tissue_list, tissue)
    }, error = function(e) {
      print(paste("Error occurred for", pheno, locus))
      print(e)
    })
  }
}

# Combine the data frames in the list into one data frame
results_df <- data.frame(pheno_locus_list, tissue_list, results_list)

# Write out the results
out_file <- paste(tissue, "coloc_ukbb_sig.txt", sep = "-")
write.table(results_df, file = file.path(out_dir, out_file), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)