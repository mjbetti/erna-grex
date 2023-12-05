library("data.table")
library("ggplot2")
library("lattice")

root_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training"

tissues <- c("Adipose_sub", "Adipose_vis", "Adrenal", "Artery_aor", "Artery_cor", "Artery_tib", "Brain_Amy", "Brain_Ant", "Brain_Cau", "Brain_Cerebellum", "Brain_Cer", "Brain_Cortex", "Brain_Fro", "Brain_Hippo", "Brain_Hypo", "Brain_Nuc", "Brain_Put", "Brain_Spi", "Brain_Sub", "Breast", "Cells_ebv", "Cells_fibroblast", "Colon_sigmoid", "Colon_transverse", "Esophagus_gast", "Esophagus_muc", "Esophagus_mus", "Heart_Atr", "Heart_Left", "Kidney", "Liver", "Lung", "Minor_salivary_gland", "Muscle_skeletal", "Nerve", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_no_sun", "Skin_sun", "Small_intestine", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_blood")

n_erna_df <- data.frame()
n_genes_df <- data.frame()

##eRNAs
for (tissue in tissues) {
	print(paste0(tissue, "..."))
	snps_array <- c()
	weights_dir <- paste0("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training/", tissue, "/weights")
	for (i in list.files(weights_dir)) {
		file <- fread(paste(weights_dir, i, sep = "/"), header = TRUE, quote = "", sep = "\t")
		df <- as.data.frame(file)
		n_snps <- nrow(df)
		snps_array <- c(snps_array, n_snps)
	}
	mean_snps <- mean(snps_array)
	new_row <- c(tissue, mean_snps)
	n_erna_df <- rbind(n_erna_df, new_row)
}

##Protein-coding genes
for (tissue in tissues) {
	print(paste0(tissue, "..."))
	n_snps_vec <- c()
	weights_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/training/tissue_snp_lists/gene_vars"
	file <- fread(paste(weights_dir, paste0(tissue, ".gene.varlist.txt"), sep = "/"), header = FALSE, sep = "|")
	df <- as.data.frame(file)
	genes <- unique(df[,2])
	for (gene in genes) {
		len <- length(which(df[,2] == gene))
		n_snps_vec <- c(n_snps_vec, len)
	}
	mean_n_snps <- mean(n_snps_vec)
	new_row <- c(tissue, mean_n_snps)
	n_genes_df <- rbind(n_genes_df, new_row)
}

n_erna_df$type <- "eRNA"
n_genes_df$type <- "Gene"

names(n_erna_df) <- c("tissue", "snps", "Category")
names(n_genes_df) <- c("tissue", "snps", "Category")

#Combine together
merged_df <- rbind(n_erna_df, n_genes_df)

pdf("erna_protein_coding_n_snps_attempt4.pdf")
ggplot(merged_df, aes(factor(tissue), as.numeric(snps), fill = Category)) + 
  geom_bar(stat="identity", position = "dodge") + 
  coord_flip() +
  scale_y_continuous(name = "Mean SNPs per gene model") +
  scale_x_discrete(name = "Tissues") +
  scale_fill_brewer(palette = "Set1")
dev.off()

print("Mean number of SNPs per eRNA model across tissues...")
print(mean(as.numeric(n_erna_df[,2])))
print("Mean number of SNPs per protein-coding gene model across tissues...")
print(mean(as.numeric(n_genes_df[,2])))

# pdf("erna_protein_coding_n_genes.pdf")
# barchart(genes~tissue,data=merged_df,groups=Category, auto.key = TRUE)
# dev.off()

# pdf("erna_protein_coding_n_genes.pdf")
# ReasonstatsDec <- Reasonstats[which(Reasonstats$Category=="Decline"),]
# ReasonstatsImp <- Reasonstats[which(Reasonstats$Category=="Improved"),]
# Reasonstats3   <- cbind(ReasonstatsImp[,3], ReasonstatsDec[,3])
# colnames(Reasonstats3) <- c("Improved", "Decline")
# rownames(Reasonstats3) <- ReasonstatsImp$Reason
# 
# windows()
#   barplot(t(Reasonstats3), beside=TRUE, ylab="number of species", 
#           cex.names=0.8, las=2, ylim=c(0,120), col=c("darkblue","red"))
#   box(bty="l")
# dev.off()
