library("data.table")

#Declare the paths of the files with the required information
gwas_path <- "/home/bettimj/aldrich_rotation/lc_gwas_saige/att5/plink/liftover_hg19/lc_afr_gwas_plink_att5.PHENO1.glm.logistic.reformat.arm.hg19.txt"
eqtl_weights_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training/Nerve/weights/ENSR00000313719.txt"
ld_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/mr/ld_scores_sccs/5.score.ld"
out_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/mr"

#Open the files as data frames
gwas_file <- fread(gwas_path, header = TRUE, sep = " ", quote = "")
gwas_df <- as.data.frame(gwas_file)

eqtl_weights_file <- fread(eqtl_weights_path, header = TRUE, sep = "\t", quote = "")
eqtl_weights_df <- as.data.frame(eqtl_weights_file)

ld_file <- fread(ld_path, header = TRUE, sep = " ", quote = "")
ld_df <- as.data.frame(ld_file)

#Merge the weights in the eQTL file with the GWAS results
merged_df <- merge(eqtl_weights_df, gwas_df, by.x = "rsid", by.y = "varid")

#Merge with LD scores
merged_df <- merge(merged_df, ld_file, by.x = "rsid", by.y = "SNP")

#Reorder and rename the relevant columns
gwas_beta <- log(merged_df$OR)
gwas_se <- sqrt(var(gwas_beta)/6998)

eqtl_beta <- merged_df$weights
eqtl_se <- sqrt(var(eqtl_beta)/328)

merged_df <- data.frame(merged_df$rsid, merged_df$effect_allele, merged_df$ldscore, eqtl_beta, eqtl_se, merged_df$p, gwas_beta, gwas_se, merged_df$P)

names(merged_df) <- c("rsid", "effect_allele_gwas", "ldscore", "eqtl_beta", "eqtl_se", "eqtl_p", "gwas_beta", "gwas_se", "gwas_p")

write.table(merged_df, file = paste(out_dir, "mrjti_ld_ENSR00000313719.txt", sep = "/"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)