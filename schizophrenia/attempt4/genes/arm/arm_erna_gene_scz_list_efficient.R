library("dplyr")

#Declare the output directory to which these results will be saved
out_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all_genes/brain_only/arm"

#Declare the path of the AALC GWAS results file
in_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all_genes/brain_only/arm/schizophrenia_eRNA_gene_brain_coord_results.tsv"

#Declare the path of the hg38 cytoBand file downloaded from UCSC. This will be used to label SNPs based on their chromosomal location (arm).
cytoband_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all_genes/brain_only/arm/cytoBand.txt.gz"

#Open the AALC concatenated output data as a data frame, and assign the p-values to a vector
in_file <- read.table(in_path, header = TRUE)
in_df <- as.data.frame(in_file)
p_vals <- in_df$pvalue

#Open the cytoBand file as a data frame
cytoband_file <- read.table(cytoband_path, header = FALSE)
cytoband_df <- as.data.frame(cytoband_file)

#Define the function that we will use to match up SNP coordinates with their cytoBand range to return chromosome arm (adapted from https://stackoverflow.com/questions/24766104/checking-if-value-in-vector-is-in-range-of-values-in-different-length-vector)
getValue <- function(x, data) {
  tmp <- data %>%
    filter(V2 <= x, x <= V3)
  return(tmp$V4)
}

#Generate a new data frame with all of the SNPs that meet our preferred p-value threshold and their corresponding p-values
cutoff_in_df <- in_df[(in_df$pvalue <= 6.327e-6),]

#Loop through and split up the coordinates by chromosome
for (chr in 1:22) {
	cutoff_in_df_chr = assign(paste0("cutoff_in_df_chr_", chr), cutoff_in_df[(cutoff_in_df$CHR == chr),])
	cytoband_df_chr = assign(paste0("cytoband_df_chr_", chr), cytoband_df[(cytoband_df[,1] == paste0("chr", chr)),])
	in_snp_coord = assign(paste0("in_snp_coord_", chr), cutoff_in_df_chr$START_POS)
	arm_vector = assign(paste0("arm_vector_chr_", chr), sapply(in_snp_coord, getValue, cytoband_df_chr))
	chr_pref = rep(chr, length(arm_vector))
	arm_vector = paste0(chr_pref, arm_vector)
	assign(paste0("cutoff_in_df_chr_", chr), cbind(cutoff_in_df_chr, arm_vector))
}

#Concatenate each of the newly arm-ammended data frames back into a single data frame
annot_cat_df <- rbind(cutoff_in_df_chr_1, cutoff_in_df_chr_2, cutoff_in_df_chr_3, cutoff_in_df_chr_4, cutoff_in_df_chr_5, cutoff_in_df_chr_6, cutoff_in_df_chr_7, cutoff_in_df_chr_8, cutoff_in_df_chr_9, cutoff_in_df_chr_10, cutoff_in_df_chr_11, cutoff_in_df_chr_12, cutoff_in_df_chr_13, cutoff_in_df_chr_14, cutoff_in_df_chr_15, cutoff_in_df_chr_16, cutoff_in_df_chr_17, cutoff_in_df_chr_18, cutoff_in_df_chr_19, cutoff_in_df_chr_20, cutoff_in_df_chr_21, cutoff_in_df_chr_22)

annot_cat_df <- annot_cat_df[order(annot_cat_df$arm_vector),]

#Slice out the columns that we want in out final output
#final_out_df <- data.frame(annot_cat_df$rsid, annot_cat_df$alternate_ids, annot_cat_df$arm, annot_cat_df$frequentist_add_pvalue)

final_out_name <- "arm_scz_erna_and_gene_twas.txt"
write.table(annot_cat_df, file = paste(out_dir, final_out_name, sep = "/"), quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)




