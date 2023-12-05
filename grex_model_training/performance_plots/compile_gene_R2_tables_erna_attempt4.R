#Declare the path of the directory containing the extracted eRNA extra tables, as well as the target output directory
tables_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training/tissue_extra_r2_lists"
out_dir <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training/tissue_extra_r2_lists/r2_tables"

#Loop through each file, extracting the gene name and R2 value, and outputting the results to a new file
for (file in list.files(tables_dir, pattern = "*.extra.predixcan.hg19.txt")) {
	tissue <- strsplit(file, "[.]")
	tissue <- unlist(sapply(tissue, "[[", 1))
	read_file <- read.table(paste(tables_dir, file, sep = "/"), header = FALSE, sep = "|", stringsAsFactors = FALSE)
	df <- as.data.frame(read_file)
	gene <- df[,1]
	r2 <- df[,3]
	new_df <- data.frame(gene, r2)
	write.table(new_df, file = paste(out_dir, paste(tissue, "r2.erna.txt", sep = "."), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}