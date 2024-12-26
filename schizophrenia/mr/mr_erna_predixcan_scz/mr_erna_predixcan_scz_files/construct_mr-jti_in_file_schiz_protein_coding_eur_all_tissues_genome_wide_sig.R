library("data.table")
library("dplyr")

#Declare the paths of the files with the required information
key_path <- "/home/bettimj/aldrich_rotation/mxcovar_keys/jti_varids_rsids_hg19_b150.txt"
gwas_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv"
eqtl_weights_root_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training"
ld_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/mr/ld_scores_gtex/all_chr.ld"
twas_results_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all_genes"
exp_dir <- "/home/bettimj/gamazon_rotation/gtexv8_exp"
out_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/mr/all_tissues/genome_wide_sig"

#Open the key, GWAS, and LD files as data frames
key_file <- fread(key_path, header = FALSE, sep = "\t", quote = "")
key_df <- as.data.frame(key_file)
key_df[,1] <- paste0(key_df[,1], "_b37")

gwas_file <- fread(gwas_path, header = TRUE, sep = "\t", quote = "")
gwas_df <- as.data.frame(gwas_file)
gwas_df$total_sample_size <- (gwas_df$NCAS + gwas_df$NCON)

ld_file <- read.table(ld_path, header = TRUE, sep = " ", stringsAsFactors = FALSE)
ld_df <- as.data.frame(ld_file)

#Loop through all of the TWAS results, concatenating them and retaining only those that are significant
#Declare the path of the hg19 cytoBand file downloaded from UCSC. This will be used to label SNPs based on their chromosomal location (arm).
cytoband_path <- "/home/bettimj/reference_genomes/cytoBand.hg19.txt.gz"

#Open the cytoBand file as a data frame
cytoband_file <- read.table(cytoband_path, header = FALSE)
cytoband_df <- as.data.frame(cytoband_file)

#Open a gene coordinates map file
map_path <- "/home/bettimj/aldrich_rotation/metaxcan/gencode.v32.GRCh37.metaxcan_map.txt"
map_file <- read.table(map_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
map_df <- as.data.frame(map_file)

#Define the function that we will use to match up SNP coordinates with their cytoBand range to return chromosome arm (adapted from https://stackoverflow.com/questions/24766104/checking-if-value-in-vector-is-in-range-of-values-in-different-length-vector)
getValue <- function(x, data) {
  tmp <- data %>%
    filter(V2 <= x, x <= V3)
  return(tmp$V4)
}

near_sig_df <- data.frame()
sig_df <- data.frame()
cat_df <- data.frame()

tissues <- c("Adipose_sub", "Adipose_vis", "Adrenal", "Artery_aor", "Artery_cor", "Artery_tib", "Brain_Amy", "Brain_Ant", "Brain_Cau", "Brain_Cerebellum", "Brain_Cer", "Brain_Cortex", "Brain_Fro", "Brain_Hippo", "Brain_Hypo", "Brain_Nuc", "Brain_Put", "Brain_Spi", "Brain_Sub", "Breast", "Cells_ebv", "Cells_fibroblast", "Colon_sigmoid", "Colon_transverse", "Esophagus_gast", "Esophagus_muc", "Esophagus_mus", "Heart_Atr", "Heart_Left", "Kidney", "Liver", "Lung", "Minor_salivary_gland", "Muscle_skeletal", "Nerve", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_no_sun", "Skin_sun", "Small_intestine", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_blood")

tissues_twas <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellum", "Brain_Cerebellar_Hemisphere", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Cultured_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")

for (i in seq(1:length(tissues))) {
    print(tissues[i])
    file <- read.csv(paste(twas_results_dir, paste0("schizophrenia2022_gene_imputations.", tissues_twas[i], ".csv"), sep = "/"), header = TRUE, stringsAsFactors = FALSE)
    df <- as.data.frame(file)
    #print(head(df))
    df$tissue <- tissues[i]
    df <- merge(df, map_df, by.x = "gene", by.y = "ENSG")
    cat_df <- rbind(cat_df, df)
    #print(head(df))
    #df <- na.omit(df[(df$pvalue <= 1.85e-6),])
    #near_sig_df <- rbind(near_sig_df, df)
    #df <- df[(df$pvalue <= 9.42e-8),]
    #sig_df <- rbind(sig_df, df)
}

n_genes <- nrow(cat_df)
sig_threshold <- 0.05/n_genes

#Loop through and split up the coordinates by chromosome
cat_df <- cat_df[(cat_df$pvalue < sig_threshold),]

for (chr in 1:22) {
    print(chr)
    in_df_sig_chr = cat_df[(cat_df$CHR == chr),]
    print(in_df_sig_chr)
    cytoband_df_chr = cytoband_df[(cytoband_df[,1] == paste0("chr", chr)),]
    in_snp_coord = cat_df$START_POS[cat_df$CHR == chr]
    arm_vector = sapply(in_snp_coord, getValue, cytoband_df_chr)
    arm_vector = arm_vector[lapply(arm_vector,length)>0]
    chr_pref = rep(paste0("chr", chr), length(arm_vector))
    arm_vector = paste0(chr_pref, arm_vector)
    print(length(arm_vector))
    print(arm_vector)
    in_df_sig_chr$locus = arm_vector
    print(in_df_sig_chr)
    assign(paste0("in_df_sig_chr_", chr), in_df_sig_chr)
}

#Concatenate each of the newly arm-ammended data frames back into a single data frame
annot_cat_sig_df <- rbind(in_df_sig_chr_1, in_df_sig_chr_2, in_df_sig_chr_3, in_df_sig_chr_4, in_df_sig_chr_5, in_df_sig_chr_6, in_df_sig_chr_7, in_df_sig_chr_8, in_df_sig_chr_9, in_df_sig_chr_10, in_df_sig_chr_11, in_df_sig_chr_12, in_df_sig_chr_13, in_df_sig_chr_14, in_df_sig_chr_15, in_df_sig_chr_16, in_df_sig_chr_17, in_df_sig_chr_18, in_df_sig_chr_19, in_df_sig_chr_20, in_df_sig_chr_21, in_df_sig_chr_22)

names(annot_cat_sig_df)[ncol(annot_cat_sig_df)] <- "arm"

final_df <- annot_cat_sig_df
final_df_ordered <- final_df[order(final_df$pvalue),]

#Loop through the expression files, compiling sample sizes for each tissue
exp_sample_sizes <- c()

for (i in seq(1:length(tissues_twas))) {
	print(paste0(tissues[i], "..."))
	file <- fread(paste(exp_dir, paste0(tissues_twas[i], ".v8.normalized_expression.bed"), sep = "/"), header = TRUE, sep = "\t", quote = "")
	df <- as.data.frame(file)
	exp_sample_size <- ncol(df) - 4
	exp_sample_sizes <- c(exp_sample_sizes, exp_sample_size)
}

exp_sample_sizes_df <- data.frame(tissues, exp_sample_sizes)

#Compile lists of the gene-tissue pairs that will be included in the MR analysis
eqtl_genes <- c()
eqtl_tissues <- c()

for (i in seq(1:nrow(final_df_ordered))) {
	eqtl_genes <- c(eqtl_genes, final_df_ordered[i,"gene"])
	eqtl_tissues <- c(eqtl_tissues, final_df_ordered[i,"tissue"])
}

#Open each of the gene files, appending their weights to a single input file (including tissue type as another column)
eqtl_cat_df <- data.frame()

for (i in seq(1:length(eqtl_genes))){
	path <- paste("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/training/tissue_snp_lists/gene_vars", paste0(eqtl_tissues[i], ".gene.varlist.txt"), sep = "/")
	file <- read.table(path, header = FALSE, sep = "|", stringsAsFactors = FALSE)
	df <- as.data.frame(file)
	df <- df[(df[,2] == eqtl_genes[i]),]
	print(eqtl_genes[i])
	print(head(df))
	df$tissue <- eqtl_tissues[i]
	df$gene <- eqtl_genes[i]
	eqtl_cat_df <- rbind(eqtl_cat_df, df)
}

#Merge the eqtl_cat_df with the respective expression model training sample sizes
eqtl_cat_df <- merge(eqtl_cat_df, exp_sample_sizes_df, by.x = "tissue", by.y = "tissues")

#Merge the eqtl_cat_df with the SNP key
merged_df <- merge(eqtl_cat_df, key_df, by.x = "V1", by.y = "V2")

#Merge the weights in the eQTL file with the GWAS results
merged_df <- merge(merged_df, gwas_df, by.x = "V1.y", by.y = "SNP")

#Merge with LD scores
merged_df <- merge(merged_df, ld_df, by.x = "V1.y", by.y = "SNP")

#Reorder and rename the relevant columns
merged_df$gwas_beta <- merged_df$BETA
merged_df$gwas_se <- sqrt(var(merged_df$gwas_beta)/merged_df$total_sample_size)

merged_df$eqtl_beta <- merged_df$V3
merged_df$eqtl_se <- sqrt(var(merged_df$eqtl_beta)/merged_df$exp_sample_sizes)

eqtl_zscore <- merged_df$eqtl_beta / merged_df$eqtl_se
merged_df$eqtl_p <- 2 * (1 - pnorm(abs(eqtl_zscore)))

merged_df <- data.frame(merged_df$V1.y, merged_df$A2, merged_df$ldscore, merged_df$eqtl_beta, merged_df$eqtl_se, merged_df$eqtl_p, merged_df$gwas_beta, merged_df$gwas_se, merged_df$PVAL, merged_df$tissue, merged_df$gene)

names(merged_df) <- c("rsid", "effect_allele_gwas", "ldscore", "eqtl_beta", "eqtl_se", "eqtl_p", "gwas_beta", "gwas_se", "gwas_p", "tissue", "gene")

print(length(unique(merged_df$gene)))
#424

write.table(merged_df, file = paste(out_dir, "mrjti_erna_all_tissues_schiz_protein_coding_2022_predixcan_genome_wide_sig.txt", sep = "/"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)