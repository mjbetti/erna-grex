library("data.table")
library("stringr")
library("dplyr")

#Declare the path of the input files
ensembl_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20161111.lifted.hg19.bed"
fantom5_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/human_permissive_enhancers_phase_1_and_2.bed.gz"
roadmap_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/roadmap/egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all_15_coreMarks_dense_enh.bed"
expids_path <- "/home/bettimj/gamazon_rotation/herna/reformatted/all.rds.erna_exp.txt"

#Open the files as data frames
ensembl_file <- fread(ensembl_path, header = FALSE, sep = "\t", quote = "")
ensembl_df <- as.data.frame(ensembl_file)
ensembl_df <- data.frame(ensembl_df[,1:5], ensembl_df[,5], ensembl_df[,5])
names(ensembl_df) <- c("chr", "left", "right", "strand", "geneid_full", "geneid", "genename")

expids_file <- fread(expids_path, header = TRUE, sep = "\t", quote = "")
expids_df <- as.data.frame(expids_file)

#A new data frame was constructed
new_df <- rbind(ensembl_df, expids_df)
new_df <- new_df[order(new_df$chr, new_df$left, new_df$right),]
new_df <- unique(new_df)
new_df <- data.frame(new_df[,1:6], "Enhancer", new_df[,7])
names(new_df)[7:8] <- c("genetype", "genename")

#Write out to a new file
out_directory <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4"
write.table(new_df, file = paste(out_directory, "ensembl87_fantom5_roadmap_hg19_reg_annos_hg19.txt", sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)