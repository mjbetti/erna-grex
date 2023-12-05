library("data.table")

#Declare the path of the original GWAS summary statistics, along with the path of the BED coordinates lifted over to hg19
gwas_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/gwas_sum_stats/schizophrenia/daner_PGC_SCZ52_0513a.hq2"
lifted_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/gwas_sum_stats/schizophrenia/liftover_hg38/schizophrenia_daner_PGC_SCZ52_0513a.hg38.bed"

#Open each of these files as data frames
gwas_file <- fread(gwas_path, header = TRUE, quote = "", sep = "\t")
gwas_df <- as.data.frame(gwas_file)
gwas_df$varid <- paste(paste0("chr", gwas_df$CHR), gwas_df$BP, gwas_df$A1, gwas_df$A2, "b37", sep = "_")

lifted_file <- fread(lifted_path, header = FALSE, quote = "", sep = "\t")
lifted_df <- as.data.frame(lifted_file)

#Merge the two files together by original variant ID (column 1 in the GWAS file and column 4 in the lifted file)
merged_df <- merge(gwas_df, lifted_df, by.x = "varid", by.y = "V4", all = TRUE)
merged_df <- data.frame(substring(merged_df$V1, 4), merged_df$V2, merged_df[,5:18])
merged_df <- data.frame(paste(paste0("chr", merged_df[,1]), merged_df[,2], merged_df$A2, merged_df$A1, "b38", sep = "_"), merged_df)
names(merged_df)[1:3] <- c("varID", "CHR", "BP")

merged_df <- merged_df[!is.na(merged_df$CHR),]

#Write the lifted summary statistics to a new file
write.table(merged_df, file = "schizophrenia_daner_PGC_SCZ52_0513a.hg38.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)