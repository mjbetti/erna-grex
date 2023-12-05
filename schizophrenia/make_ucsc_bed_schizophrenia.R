library("data.table")

#Declare the path of the AALC GWAS summary statistics, as well as the desired output directory
options(scipen = 999)
summary_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/gwas_sum_stats/schizophrenia/daner_PGC_SCZ52_0513a.hq2"
out_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/gwas_sum_stats/schizophrenia/liftover_hg38"

#Open the file as a data frame
summary_file <- fread(summary_path, quote = "", header = TRUE, sep = "\t")
summary_df <- as.data.frame(summary_file)

#Pull out the columns that will be included in the bed output
chr <- paste0("chr", summary_df$CHR)
start <- summary_df$BP
end <- (summary_df$BP + 1)
varid <- paste(chr, start, summary_df$A2, summary_df$A1, "b37", sep = "_")
new_df <- data.frame(chr, start, end, varid)
#new_df <- na.omit(new_df)

#Output as a UCSC BED file
write.table(new_df, file = paste(out_dir, "schizophrenia_daner_PGC_SCZ52_0513a.hg19.bed", sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
