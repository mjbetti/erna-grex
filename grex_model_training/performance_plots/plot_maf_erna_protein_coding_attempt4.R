library("data.table")
library("ggplot2")
library("lattice")

root_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training"

tissues <- c("Adipose_sub", "Adipose_vis", "Adrenal", "Artery_aor", "Artery_cor", "Artery_tib", "Brain_Amy", "Brain_Ant", "Brain_Cau", "Brain_Cerebellum", "Brain_Cer", "Brain_Cortex", "Brain_Fro", "Brain_Hippo", "Brain_Hypo", "Brain_Nuc", "Brain_Put", "Brain_Spi", "Brain_Sub", "Breast", "Cells_ebv", "Cells_fibroblast", "Colon_sigmoid", "Colon_transverse", "Esophagus_gast", "Esophagus_muc", "Esophagus_mus", "Heart_Atr", "Heart_Left", "Kidney", "Liver", "Lung", "Minor_salivary_gland", "Muscle_skeletal", "Nerve", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_no_sun", "Skin_sun", "Small_intestine", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_blood")

maf_path <- "/home/bettimj/gamazon_rotation/gtex8_geno/liftover_hg19/plink2.afreq"
maf_file <- fread(maf_path, header = TRUE, sep = "\t", quote = "")
maf_df <- as.data.frame(maf_file)

maf_path_rsid <- "/home/bettimj/gamazon_rotation/gtex8_geno/maf_rsid_hg19/plink2.afreq"
maf_file_rsid <- fread(maf_path_rsid, header = TRUE, sep = "\t", quote = "")
maf_df_rsid <- as.data.frame(maf_file_rsid)

#Mean
maf_erna_df <- data.frame()
maf_genes_df <- data.frame()

##eRNAs
for (tissue in tissues) {
	snps_path <- paste("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training/tissue_snp_lists", paste0(tissue, ".varlist.txt"), sep = "/")
	snps_file <- fread(snps_path, header = FALSE, sep = "|", quote = "")
	snps_df <- as.data.frame(snps_file)
	snps <- unique(snps_df[,1])
	snps_maf <- maf_df[(maf_df[,2] %in% snps),]
	mean_maf <- mean(as.numeric(snps_maf[,5]))
	new_row <- c(tissue, mean_maf)
	maf_erna_df <- rbind(maf_erna_df, new_row)
}

##Protein-coding genes
for (tissue in tissues) {
	snps_path <- paste("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/training/tissue_snp_lists/gene_vars", paste0(tissue, ".gene.varlist.txt"), sep = "/")
	snps_file <- fread(snps_path, header = FALSE, sep = "|", quote = "")
	snps_df <- as.data.frame(snps_file)
	snps <- unique(snps_df[,1])
	snps_maf <- maf_df_rsid[(maf_df_rsid[,2] %in% snps),]
	mean_maf <- mean(as.numeric(snps_maf[,5]))
	new_row <- c(tissue, mean_maf)
	maf_genes_df <- rbind(maf_genes_df, new_row)
}


maf_erna_df$type <- "eRNA"
maf_genes_df$type <- "Gene"

names(maf_erna_df) <- c("tissue", "maf", "Category")
names(maf_genes_df) <- c("tissue", "maf", "Category")

#Combine together
merged_df <- rbind(maf_erna_df, maf_genes_df)

pdf("erna_protein_coding_maf_attempt4.pdf")
ggplot(merged_df, aes(factor(tissue), as.numeric(maf), fill = Category)) + 
  geom_bar(stat="identity", position = "dodge") + 
  coord_flip() +
  scale_y_continuous(name = "MAF") +
  scale_x_discrete(name = "Tissues") +
  scale_fill_brewer(palette = "Set1")
dev.off()

#Print numbers of tissues where r2 eRNA > r2 gene and vice versa
maf_erna_df <- maf_erna_df[,1:2]
names(maf_erna_df) <- c("tissue", "maf_erna")

maf_genes_df <- maf_genes_df[,1:2]
names(maf_genes_df) <- c("tissue", "maf_genes")

merged_df <- merge(maf_erna_df, maf_genes_df, by = "tissue")

print("Number tissues eRNA MAF greater than gene MAF:")
print(nrow(merged_df[(merged_df$maf_erna > merged_df$maf_genes),]))

print("Number tissues eRNA MAF less than gene MAF:")
print(nrow(merged_df[(merged_df$maf_erna < merged_df$maf_genes),]))

#Median
maf_erna_df <- data.frame()
maf_genes_df <- data.frame()

##eRNAs
for (tissue in tissues) {
	snps_path <- paste("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training/tissue_snp_lists", paste0(tissue, ".varlist.txt"), sep = "/")
	snps_file <- fread(snps_path, header = FALSE, sep = "|", quote = "")
	snps_df <- as.data.frame(snps_file)
	snps <- unique(snps_df[,1])
	snps_maf <- maf_df[(maf_df[,2] %in% snps),]
	mean_maf <- median(as.numeric(snps_maf[,5]))
	new_row <- c(tissue, mean_maf)
	maf_erna_df <- rbind(maf_erna_df, new_row)
}

##Protein-coding genes
for (tissue in tissues) {
	snps_path <- paste("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/training/tissue_snp_lists/gene_vars", paste0(tissue, ".gene.varlist.txt"), sep = "/")
	snps_file <- fread(snps_path, header = FALSE, sep = "|", quote = "")
	snps_df <- as.data.frame(snps_file)
	snps <- unique(snps_df[,1])
	snps_maf <- maf_df_rsid[(maf_df_rsid[,2] %in% snps),]
	mean_maf <- median(as.numeric(snps_maf[,5]))
	new_row <- c(tissue, mean_maf)
	maf_genes_df <- rbind(maf_genes_df, new_row)
}


maf_erna_df$type <- "eRNA"
maf_genes_df$type <- "Gene"

names(maf_erna_df) <- c("tissue", "maf", "Category")
names(maf_genes_df) <- c("tissue", "maf", "Category")

#Combine together
merged_df <- rbind(maf_erna_df, maf_genes_df)

pdf("erna_protein_coding_median_maf_attempt4.pdf")
ggplot(merged_df, aes(factor(tissue), as.numeric(maf), fill = Category)) + 
  geom_bar(stat="identity", position = "dodge") + 
  coord_flip() +
  scale_y_continuous(name = "MAF") +
  scale_x_discrete(name = "Tissues") +
  scale_fill_brewer(palette = "Set1")
dev.off()

#Print numbers of tissues where r2 eRNA > r2 gene and vice versa
maf_erna_df <- maf_erna_df[,1:2]
names(maf_erna_df) <- c("tissue", "maf_erna")

maf_genes_df <- maf_genes_df[,1:2]
names(maf_genes_df) <- c("tissue", "maf_genes")

merged_df <- merge(maf_erna_df, maf_genes_df, by = "tissue")

print("Number tissues eRNA MAF greater than gene MAF:")
print(nrow(merged_df[(merged_df$maf_erna > merged_df$maf_genes),]))

print("Number tissues eRNA MAF less than gene MAF:")
print(nrow(merged_df[(merged_df$maf_erna < merged_df$maf_genes),]))


#Compare MAFs in each tissue using Wilcoxon test
for (tissue in tissues) {
	print(paste(tissue, "..."))
	snps_path_erna <- paste("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training/tissue_snp_lists", paste0(tissue, ".varlist.txt"), sep = "/")
	snps_file_erna <- fread(snps_path_erna, header = FALSE, sep = "|", quote = "")
	snps_df_erna <- as.data.frame(snps_file_erna)
	snps_erna <- unique(snps_df_erna[,1])
	snps_maf_erna <- maf_df[(maf_df[,2] %in% snps_erna),]
	
	snps_path_genes <- paste("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/training/tissue_snp_lists/gene_vars", paste0(tissue, ".gene.varlist.txt"), sep = "/")
	snps_file_genes <- fread(snps_path_genes, header = FALSE, sep = "|", quote = "")
	snps_df_genes <- as.data.frame(snps_file_genes)
	snps_genes <- unique(snps_df_genes[,1])
	snps_maf_genes <- maf_df_rsid[(maf_df_rsid[,2] %in% snps_genes),]
	
	print(wilcox.test(snps_maf_erna[,5], snps_maf_genes[,5]))
}