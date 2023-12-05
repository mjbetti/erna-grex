library("data.table")

sig_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all_genes/brain_only/arm/arm_scz_erna_and_gene_twas.txt"

#Open the associations file as a data frame
sig_file <- fread(sig_path, header = TRUE, sep = " ", quote = "")
sig_df <- as.data.frame(sig_file)

for (line in seq(1, nrow(sig_df))) {
	if (sig_df[line,"tissue"] == "Adipose_sub") {
 		sig_df[line,"tissue"] <- "Adipose_Subcutaneous"
	} else if (sig_df[line,"tissue"] == "Adipose_vis") {
		sig_df[line,"tissue"] <- "Adipose_Visceral_Omentum"
	} else if (sig_df[line,"tissue"] == "Adrenal") {
		sig_df[line,"tissue"] <- "Adrenal_Gland"
	} else if (sig_df[line,"tissue"] == "Artery_aor") {
		sig_df[line,"tissue"] <- "Artery_Aorta"
	} else if (sig_df[line,"tissue"] == "Artery_cor") {
		sig_df[line,"tissue"] <- "Artery_Coronary"
	} else if (sig_df[line,"tissue"] == "Artery_tib") {
		sig_df[line,"tissue"] <- "Artery_Tibial"
	} else if (sig_df[line,"tissue"] == "Brain_Amy") {
		sig_df[line,"tissue"] <- "Brain_Amygdala"
	} else if (sig_df[line,"tissue"] == "Brain_Ant") {
		sig_df[line,"tissue"] <- "Brain_Anterior_cingulate_cortex_BA24"
	} else if (sig_df[line,"tissue"] == "Brain_Cau") {
		sig_df[line,"tissue"] <- "Brain_Caudate_basal_ganglia"
	} else if (sig_df[line,"tissue"] == "Brain_Cerebellum") {
		sig_df[line,"tissue"] <- "Brain_Cerebellum"
	} else if (sig_df[line,"tissue"] == "Brain_Cer") {
		sig_df[line,"tissue"] <- "Brain_Cerebellar_Hemisphere"
	} else if (sig_df[line,"tissue"] == "Brain_Cortex") {
		sig_df[line,"tissue"] <- "Brain_Cortex"
	} else if (sig_df[line,"tissue"] == "Brain_Fro") {
		sig_df[line,"tissue"] <- "Brain_Frontal_Cortex_BA9"
	} else if (sig_df[line,"tissue"] == "Brain_Hippo") {
		sig_df[line,"tissue"] <- "Brain_Hippocampus"
	} else if (sig_df[line,"tissue"] == "Brain_Hypo") {
		sig_df[line,"tissue"] <- "Brain_Hypothalamus"
	} else if (sig_df[line,"tissue"] == "Brain_Nuc") {
		sig_df[line,"tissue"] <- "Brain_Nucleus_accumbens_basal_ganglia"
	} else if (sig_df[line,"tissue"] == "Brain_Put") {
		sig_df[line,"tissue"] <- "Brain_Putamen_basal_ganglia"
	} else if (sig_df[line,"tissue"] == "Brain_Spi") {
		sig_df[line,"tissue"] <- "Brain_Spinal_cord_cervical_c-1"
	} else if (sig_df[line,"tissue"] == "Brain_Sub") {
		sig_df[line,"tissue"] <- "Brain_Substantia_nigra"
	} else if (sig_df[line,"tissue"] == "Breast") {
		sig_df[line,"tissue"] <- "Breast_Mammary_Tissue"
	} else if (sig_df[line,"tissue"] == "Cells_ebv") {
		sig_df[line,"tissue"] <- "Cells_EBV-transformed_lymphocytes"
	} else if (sig_df[line,"tissue"] == "Cells_fibroblast") {
		sig_df[line,"tissue"] <- "Cells_Cultured_fibroblasts"
	} else if (sig_df[line,"tissue"] == "Colon_sigmoid") {
		sig_df[line,"tissue"] <- "Colon_Sigmoid"
	} else if (sig_df[line,"tissue"] == "Colon_transverse") {
		sig_df[line,"tissue"] <- "Colon_Transverse"
	} else if (sig_df[line,"tissue"] == "Esophagus_gast") {
		sig_df[line,"tissue"] <- "Esophagus_Gastroesophageal_Junction"
	} else if (sig_df[line,"tissue"] == "Esophagus_muc") {
		sig_df[line,"tissue"] <- "Esophagus_Mucosa"
	} else if (sig_df[line,"tissue"] == "Esophagus_mus") {
		sig_df[line,"tissue"] <- "Esophagus_Muscularis"
	} else if (sig_df[line,"tissue"] == "Heart_Atr") {
		sig_df[line,"tissue"] <- "Heart_Atrial_Appendage"
	} else if (sig_df[line,"tissue"] == "Heart_Left") {
		sig_df[line,"tissue"] <- "Heart_Left_Ventricle"
	} else if (sig_df[line,"tissue"] == "Kidney") {
		sig_df[line,"tissue"] <- "Kidney_Cortex"
	} else if (sig_df[line,"tissue"] == "Liver") {
		sig_df[line,"tissue"] <- "Liver"
	} else if (sig_df[line,"tissue"] == "Lung") {
		sig_df[line,"tissue"] <- "Lung"
	} else if (sig_df[line,"tissue"] == "Minor_salivary_gland") {
		sig_df[line,"tissue"] <- "Minor_Salivary_Gland"
	} else if (sig_df[line,"tissue"] == "Muscle_skeletal") {
		sig_df[line,"tissue"] <- "Muscle_Skeletal"
	} else if (sig_df[line,"tissue"] == "Nerve") {
		sig_df[line,"tissue"] <- "Nerve_Tibial"
	} else if (sig_df[line,"tissue"] == "Ovary") {
		sig_df[line,"tissue"] <- "Ovary"
	} else if (sig_df[line,"tissue"] == "Pancreas") {
		sig_df[line,"tissue"] <- "Pancreas"
	} else if (sig_df[line,"tissue"] == "Pituitary") {
		sig_df[line,"tissue"] <- "Pituitary"
	} else if (sig_df[line,"tissue"] == "Prostate") {
		sig_df[line,"tissue"] <- "Prostate"
	} else if (sig_df[line,"tissue"] == "Skin_no_sun") {
		sig_df[line,"tissue"] <- "Skin_Not_Sun_Exposed_Suprapubic"
	} else if (sig_df[line,"tissue"] == "Skin_sun") {
		sig_df[line,"tissue"] <- "Skin_Sun_Exposed_Lower_leg"
	} else if (sig_df[line,"tissue"] == "Small_intestine") {
		sig_df[line,"tissue"] <- "Small_Intestine_Terminal_Ileum"
	} else if (sig_df[line,"tissue"] == "Spleen") {
		sig_df[line,"tissue"] <- "Spleen"
	} else if (sig_df[line,"tissue"] == "Stomach") {
		sig_df[line,"tissue"] <- "Stomach"
	} else if (sig_df[line,"tissue"] == "Testis") {
		sig_df[line,"tissue"] <- "Testis"
	} else if (sig_df[line,"tissue"] == "Thyroid") {
		sig_df[line,"tissue"] <- "Thyroid"
	} else if (sig_df[line,"tissue"] == "Uterus") {
		sig_df[line,"tissue"] <- "Uterus"
	} else if (sig_df[line,"tissue"] == "Vagina") {
		sig_df[line,"tissue"] <- "Vagina"
	} else if (sig_df[line,"tissue"] == "Whole_blood") {
		sig_df[line,"tissue"] <- "Whole_Blood"
	}
}

#Remove those without at least one eRNA and gene association at the same locus in the same brain tissue
#sig_df <- na.omit(sig_df[duplicated(c(sig_df$tissue, sig_df$arm_vector)),])

erna_df <- sig_df[startsWith(sig_df$gene, "ENSR"),]
erna_df$tissuearm <- paste(erna_df$tissue, erna_df$arm_vector, sep = "_")

gene_df <- sig_df[startsWith(sig_df$gene, "ENSG"),]
gene_df$tissuearm <- paste(gene_df$tissue, gene_df$arm_vector, sep = "_")

erna_df <- erna_df[(erna_df$tissuearm %in% gene_df$tissuearm),]
gene_df <- gene_df[(gene_df$tissuearm %in% erna_df$tissuearm),]

#Compile a new data frame with all loci with at least one eRNA and one gene association per locus
all_df <- rbind(erna_df, gene_df)
all_df <- all_df[order(all_df$arm_vector, all_df$tissue),]

all_df <- all_df[duplicated(all_df$tissuearm),] 

erna_df <- all_df[startsWith(all_df$gene, "ENSR"),]
gene_df <- all_df[startsWith(all_df$gene, "ENSG"),]

erna_df <- erna_df[(erna_df$tissuearm %in% gene_df$tissuearm),]
gene_df <- gene_df[(gene_df$tissuearm %in% erna_df$tissuearm),]

all_df <- rbind(erna_df, gene_df)
all_df <- all_df[order(all_df$tissuearm),]

tissuearm <- all_df$tissuearm
tab <- table(tissuearm)
tissuearm <- tissuearm[tissuearm %in% names(tab[tab == 2])]

pair_df <- all_df[(all_df$tissuearm %in% tissuearm),]

#For the paired associations, find the Spearman correlation of the expression values between the two
brain_cerebellum_erna_path <- "/home/bettimj/gamazon_rotation/herna/normalized_to_covs/Brain_Cerebellum.rds.erna_exp.norm.sex_age_platform_pcs_peer1_10.txt"
brain_cerebellum_gene_path <- "/home/bettimj/gamazon_rotation/gtexv8_exp/Brain_Cerebellum.v8.normalized_expression.bed"

brain_fro_erna_path <- "/home/bettimj/gamazon_rotation/herna/normalized_to_covs/Brain_Fro.rds.erna_exp.norm.sex_age_platform_pcs_peer1_10.txt"
brain_fro_gene_path <- "/home/bettimj/gamazon_rotation/gtexv8_exp/Brain_Frontal_Cortex_BA9.v8.normalized_expression.bed"

brain_cerebellum_erna_file <- read.table(brain_cerebellum_erna_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
brain_cerebellum_erna_df <- as.data.frame(brain_cerebellum_erna_file)

brain_cerebellum_gene_file <- read.table(brain_cerebellum_gene_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
brain_cerebellum_gene_df <- as.data.frame(brain_cerebellum_gene_file)

brain_fro_erna_file <- read.table(brain_fro_erna_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
brain_fro_erna_df <- as.data.frame(brain_fro_erna_file)

brain_fro_gene_file <- read.table(brain_fro_gene_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
brain_fro_gene_df <- as.data.frame(brain_fro_gene_file)

erna_df_cer <- brain_cerebellum_erna_df[(brain_cerebellum_erna_df$geneid == "ENSR00000116923"),]
sample_names <- names(erna_df_cer[,2:ncol(erna_df_cer)])
erna_cer_values <- as.numeric(erna_df_cer[1,2:ncol(erna_df_cer)])
new_erna_df_cer <- data.frame(sample_names, erna_cer_values)

erna_df_fro <- brain_fro_erna_df[(brain_fro_erna_df$geneid == "ENSR00000128482"),]
sample_names <- names(brain_fro_erna_df[,2:ncol(brain_fro_erna_df)])
erna_fro_values <- as.numeric(erna_df_fro[1,2:ncol(erna_df_fro)])
new_erna_df_fro <- data.frame(sample_names, erna_fro_values)

gene_df_cer <- brain_cerebellum_gene_df[startsWith(brain_cerebellum_gene_df$gene_id, "ENSG00000170802"),]
sample_names <- names(gene_df_cer[,5:ncol(gene_df_cer)])
gene_cer_values <- as.numeric(gene_df_cer[1,5:ncol(gene_df_cer)])
new_gene_df_cer <- data.frame(sample_names, gene_cer_values)

gene_df_fro <- brain_fro_gene_df[startsWith(brain_fro_gene_df$gene_id, "ENSG00000162971"),]
sample_names <- names(brain_fro_gene_df[,5:ncol(brain_fro_gene_df)])
gene_fro_values <- as.numeric(gene_df_fro[1,5:ncol(gene_df_fro)])
new_gene_df_fro <- data.frame(sample_names, gene_fro_values)

merged_df_cer <- merge(new_erna_df_cer, new_gene_df_cer, by = "sample_names")
merged_df_fro <- merge(new_erna_df_fro, new_gene_df_fro, by = "sample_names")

merged_df_cer.cor = cor(merged_df_cer[,2:ncol(merged_df_cer)], method = c("spearman"))
merged_df_cer_corrval <- merged_df_cer.cor[1,2]

merged_df_fro.cor = cor(merged_df_fro[,2:ncol(merged_df_fro)], method = c("spearman"))
merged_df_fro_corrval <- merged_df_fro.cor[1,2]

pair_df$corr[pair_df$tissue == "Brain_Cerebellum"] <- merged_df_cer_corrval
pair_df$corr[pair_df$tissue == "Brain_Frontal_Cortex_BA9"] <- merged_df_fro_corrval

#Write out the all dataframe and pair dataframes
out_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all_genes/brain_only/arm"

write.table(all_df, file = paste(out_dir, "erna_gene_association_loci_arm_scz_erna_and_gene_twas.txt", sep = "/"), quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)

write.table(pair_df, file = paste(out_dir, "paired_erna_gene_association_loci_arm_scz_erna_and_gene_twas.txt", sep = "/"), quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)

