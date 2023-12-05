library("data.table")
library("stringr")
library("dplyr")

#Declare the path of the input GFF file
in_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20161117.gff.gz"

#Open the file as a data frame
in_file <- fread(in_path, header = FALSE, sep = "\t", quote = "")
in_df <- as.data.frame(in_file)

#Split up the last column to extract the ENSR IDs
id_strs_full <- in_df[,9]
id_strs_full <- strsplit(id_strs_full, ";")
id_strs_full <- unlist(lapply(id_strs_full, `[[`, 1))
id_strs_full <- strsplit(as.character(id_strs_full), "=")
ensr_ids <- unlist(lapply(id_strs_full, `[[`, 2))

#A new data frame was constructed
new_df <- data.frame(paste0("chr", in_df[,1]), in_df[,4], in_df[,5], in_df[,6], ensr_ids, ensr_ids, ensr_ids)
names(new_df) <- c("chr", "left", "right", "strand", "geneid", "genetype", "genename")

#Write out to a new file
out_directory <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4"
write.table(new_df, file = paste(out_directory, "homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20161117.bed", sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)