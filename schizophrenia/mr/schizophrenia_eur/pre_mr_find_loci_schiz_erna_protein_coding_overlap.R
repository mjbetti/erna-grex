library("data.table")
library("dplyr")

##eRNA
#Declare the paths of the files with the required information
gwas_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv"
eqtl_weights_root_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training"
ld_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/mr/ld_scores_gtex/all_chr.ld"
twas_results_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all"
exp_dir <- "/home/bettimj/gamazon_rotation/herna/reformatted"
out_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/mr"

#Open the GWAS and LD files as data frames
gwas_file <- fread(gwas_path, header = TRUE, sep = "\t", quote = "")
gwas_df <- as.data.frame(gwas_file)
gwas_df$total_sample_size <- (gwas_df$NCAS + gwas_df$NCON)

ld_file <- read.table(ld_path, header = TRUE, sep = " ", stringsAsFactors = FALSE)
ld_df <- as.data.frame(ld_file)

#Loop through all of the TWAS results, concatenating them and retaining only those that are significant
#Declare the path of the hg38 cytoBand file downloaded from UCSC. This will be used to label SNPs based on their chromosomal location (arm).
cytoband_path <- "/home/bettimj/aldrich_rotation/iges2022/manhattan/cytoBand.txt.gz"

#Open the cytoBand file as a data frame
cytoband_file <- read.table(cytoband_path, header = FALSE)
cytoband_df <- as.data.frame(cytoband_file)

#Open a gene coordinates map file
map_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/ensembl87_fantom5_roadmap_hg19_reg_annos_hg19_map.txt"
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

tissues <- c("Brain_Amy", "Brain_Ant", "Brain_Cau", "Brain_Cerebellum", "Brain_Cer", "Brain_Cortex", "Brain_Fro", "Brain_Hippo", "Brain_Hypo", "Brain_Nuc", "Brain_Put", "Brain_Spi", "Brain_Sub")

for (tissue in tissues) {
    print(tissue)
    file <- read.csv(paste(twas_results_dir, paste0("schizophrenia2022_eRNA_imputations.", tissue, ".csv"), sep = "/"), header = TRUE, stringsAsFactors = FALSE)
    df <- as.data.frame(file)
    #print(head(df))
    df$tissue <- tissue
    df <- merge(df, map_df, by.x = "gene", by.y = "ENSG")
    cat_df <- rbind(cat_df, df)
    #print(head(df))
    #df <- na.omit(df[(df$pvalue <= 1.85e-6),])
    #near_sig_df <- rbind(near_sig_df, df)
    #df <- df[(df$pvalue <= 9.42e-8),]
    #sig_df <- rbind(sig_df, df)
}

n_genes <- length(unique(cat_df$gene))
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
final_df_ordered_erna <- final_df[order(final_df$pvalue),]

##Protein-coding genes
key_path <- "/home/bettimj/aldrich_rotation/mxcovar_keys/jti_varids_rsids_hg19_b150.txt"
gwas_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv"
eqtl_weights_root_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training"
ld_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/mr/ld_scores_gtex/all_chr.ld"
twas_results_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all_genes"
exp_dir <- "/home/bettimj/gamazon_rotation/gtexv8_exp"
out_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/mr"

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
#Declare the path of the hg38 cytoBand file downloaded from UCSC. This will be used to label SNPs based on their chromosomal location (arm).
cytoband_path <- "/home/bettimj/aldrich_rotation/iges2022/manhattan/cytoBand.txt.gz"

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

tissues <- c("Brain_Amy", "Brain_Ant", "Brain_Cau", "Brain_Cerebellum", "Brain_Cer", "Brain_Cortex", "Brain_Fro", "Brain_Hippo", "Brain_Hypo", "Brain_Nuc", "Brain_Put", "Brain_Spi", "Brain_Sub")

tissues_twas <- c("Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellum", "Brain_Cerebellar_Hemisphere", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra")

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

n_genes <- length(unique(cat_df$gene))
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
final_df_ordered_protein_coding <- final_df[order(final_df$pvalue),]