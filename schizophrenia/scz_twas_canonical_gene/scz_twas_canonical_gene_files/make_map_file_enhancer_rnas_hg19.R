library("data.table")

#Declare the path of the GENCODE hg38 gene annotation file
anno_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/ensembl87_fantom5_roadmap_hg19_reg_annos_hg19.txt"

#Open the file as a data frame
anno_file <- fread(anno_path, header = TRUE, sep = "\t", quote = "")
anno_df <- as.data.frame(anno_file)

subchr <- substring(anno_df$chr, 4)

#Reformat to make a proper map file
new_df <- data.frame(anno_df$genename, anno_df$geneid, subchr, anno_df$left, anno_df$right)
names(new_df) <- c("Gene", "ENSG", "CHR", "START_POS", "END_POS")

#Write out as a new map file
write.table(new_df, file = "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/ensembl87_fantom5_roadmap_hg19_reg_annos_hg19_map.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)