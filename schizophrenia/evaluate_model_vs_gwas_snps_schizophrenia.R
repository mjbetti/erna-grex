library("data.table")

#Declare the paths of the lifted over schizophrenia GWAS summary statistics and the SNPs in one of the brain eRNA models. I chose cerebellum since it has the largest number of SNPs
gwas_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/gwas_sum_stats/schizophrenia/liftover_hg38/schizophrenia_daner_PGC_SCZ52_0513a.hg38.txt"
model_snps_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/tissue_snp_lists/Brain_Cerebellum.varlist.txt"

#Open each of the files as data frames
gwas_file <- fread(gwas_path, header = TRUE, quote = "")
gwas_df <- as.data.frame(gwas_file)

model_snps_file <- fread(model_snps_path, header = FALSE, sep = "|", quote = "")
model_snps_df <- as.data.frame(model_snps_file)

model_in_gwas <- model_snps_df[(model_snps_df$V1 %in% gwas_df$varID),]
print(nrow(model_in_gwas))