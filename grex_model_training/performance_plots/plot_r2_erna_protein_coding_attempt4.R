library("data.table")
library("ggplot2")
library("lattice")

root_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training"

tissues <- c("Adipose_sub", "Adipose_vis", "Adrenal", "Artery_aor", "Artery_cor", "Artery_tib", "Brain_Amy", "Brain_Ant", "Brain_Cau", "Brain_Cerebellum", "Brain_Cer", "Brain_Cortex", "Brain_Fro", "Brain_Hippo", "Brain_Hypo", "Brain_Nuc", "Brain_Put", "Brain_Spi", "Brain_Sub", "Breast", "Cells_ebv", "Cells_fibroblast", "Colon_sigmoid", "Colon_transverse", "Esophagus_gast", "Esophagus_muc", "Esophagus_mus", "Heart_Atr", "Heart_Left", "Kidney", "Liver", "Lung", "Minor_salivary_gland", "Muscle_skeletal", "Nerve", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_no_sun", "Skin_sun", "Small_intestine", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_blood")

r2_erna_df <- data.frame()
r2_genes_df <- data.frame()

##eRNAs
for (tissue in tissues) {
	r2_path <- paste("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training/tissue_extra_r2_lists/r2_tables", paste0(tissue, ".r2.erna.txt"), sep = "/")
	r2_file <- fread(r2_path, header = FALSE, sep = "\t", quote = "")
	r2_df <- as.data.frame(r2_file)
	mean_r2 <- mean(r2_df[,2])
	new_row <- c(tissue, mean_r2)
	r2_erna_df <- rbind(r2_erna_df, new_row)
}

##Protein-coding genes
for (tissue in tissues) {
	r2_path <- paste("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/training/tissue_extra_r2_lists/gene_vars/r2_tables", paste0(tissue, ".r2.gene.txt"), sep = "/")
	r2_file <- fread(r2_path, header = FALSE, sep = "\t", quote = "")
	r2_df <- as.data.frame(r2_file)
	mean_r2 <- mean(r2_df[,2])
	new_row <- c(tissue, mean_r2)
	r2_genes_df <- rbind(r2_genes_df, new_row)
}

r2_erna_df$type <- "eRNA"
r2_genes_df$type <- "Gene"

names(r2_erna_df) <- c("tissue", "r2", "Category")
names(r2_genes_df) <- c("tissue", "r2", "Category")

#Combine together
merged_df <- rbind(r2_erna_df, r2_genes_df)

pdf("erna_protein_coding_r2_attempt4.pdf")
ggplot(merged_df, aes(factor(tissue), as.numeric(r2), fill = Category)) + 
  geom_bar(stat="identity", position = "dodge") + 
  coord_flip() +
  scale_y_continuous(name = "Predictive r2") +
  scale_x_discrete(name = "Tissues") +
  scale_fill_brewer(palette = "Set1")
dev.off()

#Print numbers of tissues where r2 eRNA > r2 gene and vice versa
r2_erna_df <- r2_erna_df[,1:2]
names(r2_erna_df) <- c("tissue", "r2_erna")

r2_genes_df <- r2_genes_df[,1:2]
names(r2_genes_df) <- c("tissue", "r2_genes")

merged_df <- merge(r2_erna_df, r2_genes_df, by = "tissue")

print("Number tissues eRNA r2 greater than gene r2:")
print(nrow(merged_df[(merged_df$r2_erna > merged_df$r2_genes),]))

print("Number tissues eRNA r2 less than gene r2:")
print(nrow(merged_df[(merged_df$r2_erna < merged_df$r2_genes),]))


#Compare r2s in each tissue using Wilcoxon test
for (tissue in tissues) {
	print(paste0(tissue, "..."))
	r2_path_erna <- paste("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training/tissue_extra_r2_lists/r2_tables", paste0(tissue, ".r2.erna.txt"), sep = "/")
	r2_file_erna <- fread(r2_path_erna, header = FALSE, sep = "\t", quote = "")
	r2_df_erna <- as.data.frame(r2_file_erna)
	
	r2_path_gene <- paste("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/training/tissue_extra_r2_lists/gene_vars/r2_tables", paste0(tissue, ".r2.gene.txt"), sep = "/")
	r2_file_gene <- fread(r2_path_gene, header = FALSE, sep = "\t", quote = "")
	r2_df_gene <- as.data.frame(r2_file_gene)
	
	print(wilcox.test(r2_df_erna[,2], r2_df_gene[,2]))
}