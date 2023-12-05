library("data.table")

#Declare the path of the root directory with associations, as well as the output directory
root_dir <- "/data/g_gamazon_lab/UKBioBank/SprediXcan/eRNA_Results"
out_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/uk_biobank"

all_df_sig <- data.frame()
all_df_near_sig <- data.frame()

for (tissue in list.files(root_dir)) {
	for (twas in list.files(paste(root_dir, tissue, sep = "/"))) {
		print(paste(tissue, twas))
		file <- fread(paste(root_dir, tissue, twas, sep = "/"), header = TRUE, sep = ",", quote = "")
		df <- as.data.frame(file)
		df_sig <- df[df$pvalue <= 1.2169e-6,]
		df_near_sig <- df[df$pvalue <= 3.4437e-6,]
		tissue_pheno <- strsplit(twas, "[.]")
		tissue_pheno <- unlist(sapply(tissue_pheno, "[[", 1))
		if (nrow(df_sig) > 0) {
		df_sig$tissue_pheno <- tissue_pheno
		all_df_sig <- rbind(all_df_sig, df_sig)
		}
		if (nrow(df_near_sig) > 0) {
		df_near_sig$tissue_pheno <- tissue_pheno
		all_df_near_sig <- rbind(all_df_near_sig, df_near_sig)
		}
	}
}

all_df_sig[1:4574,13] <- paste0("Adipose_sub_", all_df_sig[1:4574,13])

all_df_sig_ordered <- all_df_sig[order(all_df_sig$pvalue),]
all_df_near_sig_ordered <- all_df_near_sig[order(all_df_near_sig$pvalue),]

#Write out the ordered data frames as new files
write.table(all_df_sig_ordered, file = paste(out_dir, "ukbb_erna_all_significant.csv", sep = "/"), quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

write.table(all_df_near_sig_ordered, file = paste(out_dir, "ukbb_erna_all_near_significant.csv", sep = "/"), quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

#For the genome-wide significant eRNAs, print the tissue and phenotype with the most significant associations
##Tissue
tissues <- c("Adipose_sub", "Adipose_vis", "Adrenal", "Artery_aor", "Artery_cor", "Artery_tib", "Brain_Amy", "Brain_Ant", "Brain_Cau", "Brain_Cerebellum", "Brain_Cer", "Brain_Cortex", "Brain_Fro", "Brain_Hippo", "Brain_Hypo", "Brain_Nuc", "Brain_Put", "Brain_Spi", "Brain_Sub", "Breast", "Cells_ebv", "Cells_fibroblast", "Colon_sigmoid", "Colon_transverse", "Esophagus_gast", "Esophagus_muc", "Esophagus_mus", "Heart_Atr", "Heart_Left", "Kidney", "Liver", "Lung", "Minor_salivary_gland", "Muscle_skeletal", "Nerve", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_no_sun", "Skin_sun", "Small_intestine", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_blood")

tissue_counts <- c()
for (tissue in tissues) {
	tissue_counts <- c(tissue_counts, nrow(all_df_sig[startsWith(all_df_sig[,13], tissue),]))
}

tissue_df <- data.frame(tissues, tissue_counts)
write.table(tissue_df, file = paste(out_dir, "ukbb_tissue_counts_genome_sig.txt", sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

new_names <- c()
##Phenotype
for (tissue in tissues) {
	name <- all_df_sig[startsWith(all_df_sig[,13], tissue),13]
	name_no_tissue <- (substring(name, nchar(tissue) + 1))
	new_names <- name_no_tissue
	#name_no_tissue[startsWith(name_no_tissue, "ebellum")] <- (substring(name, 7))
	#name_no_tissue[startsWith(name_no_tissue, "_")] <- (substring(name, 1))
}

new_names[startsWith(new_names, "ebellum")] <- (substring(new_names, 7))
new_names[startsWith(new_names, "_")] <- (substring(new_names, 1))

pheno_df <- as.data.frame(table(new_names))
pheno_df <- pheno_df[order(pheno_df[,2]),]

write.table(pheno_df, file = paste(out_dir, "ukbb_pheno_association_counts_genome_sig.txt", sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#For the suggested significant eRNAs, print the tissue and phenotype with the most significant associations
##Tissue
tissues <- c("Adipose_sub", "Adipose_vis", "Adrenal", "Artery_aor", "Artery_cor", "Artery_tib", "Brain_Amy", "Brain_Ant", "Brain_Cau", "Brain_Cerebellum", "Brain_Cer", "Brain_Cortex", "Brain_Fro", "Brain_Hippo", "Brain_Hypo", "Brain_Nuc", "Brain_Put", "Brain_Spi", "Brain_Sub", "Breast", "Cells_ebv", "Cells_fibroblast", "Colon_sigmoid", "Colon_transverse", "Esophagus_gast", "Esophagus_muc", "Esophagus_mus", "Heart_Atr", "Heart_Left", "Kidney", "Liver", "Lung", "Minor_salivary_gland", "Muscle_skeletal", "Nerve", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_no_sun", "Skin_sun", "Small_intestine", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_blood")

tissue_counts <- c()
for (tissue in tissues) {
	tissue_counts <- c(tissue_counts, nrow(all_df_near_sig[startsWith(all_df_near_sig[,13], tissue),]))
}

tissue_df <- data.frame(tissues, tissue_counts)
write.table(tissue_df, file = paste(out_dir, "ukbb_tissue_counts_genome_near_sig.txt", sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

new_names <- c()
##Phenotype
for (tissue in tissues) {
	name <- all_df_near_sig[startsWith(all_df_near_sig[,13], tissue),13]
	name_no_tissue <- (substring(name, nchar(tissue) + 1))
	new_names <- name_no_tissue
	#name_no_tissue[startsWith(name_no_tissue, "ebellum")] <- (substring(name, 7))
	#name_no_tissue[startsWith(name_no_tissue, "_")] <- (substring(name, 1))
}

new_names[startsWith(new_names, "ebellum")] <- (substring(new_names, 7))
new_names[startsWith(new_names, "_")] <- (substring(new_names, 1))

pheno_df <- as.data.frame(table(new_names))
pheno_df <- pheno_df[order(pheno_df[,2]),]

write.table(pheno_df, file = paste(out_dir, "ukbb_pheno_association_counts_genome_near_sig.txt", sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)