library("data.table")
library("optparse")

option_list = list(
	make_option(c("-i", "--in_file_name"), type = "character", default = NULL, help = "name of the input file", metavar = "character"),
	make_option(c("-d", "--in_dir"), type = "character", default = NULL, help = "directory of the input file", metavar = "character"),
	make_option(c("-l", "--loci_path"), type = "character", default = NULL, help = "path of file containing chromosomal loci", metavar = "character")
	#make_option(c("-o", "--output_path"), type = "character", default = NULL, help = "path of the output directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

in_file_name <- opt$in_file_name
in_dir <- opt$in_dir
loci_path <- opt$loci_path
#out_path <- opt$out_path
out_path <- "/home/bettimj/gamazon_rotation/ukbb_with_arm/coloc_gwas_in"

#in_dir <- "/data/g_gamazon_lab/UKBioBank/both_sexes"
#out_dir <- "/home/bettimj/gamazon_rotation/ukbb_with_arm"

pheno <- strsplit(in_file_name, "[.]")
pheno <- unlist(lapply(pheno, `[`, 1))

input_file <- fread(paste(in_dir, in_file_name, sep = "/"), quote = FALSE, sep = "\t", header = TRUE)
input_df <- as.data.frame(input_file)

loci_file <- fread(loci_path, quote = FALSE, sep = "\t", header = TRUE)
loci_df <- as.data.frame(loci_file)

merged_df <- merge(input_df, loci_df, by = "variant")

sig_df <- merged_df[(merged_df$pval < 5e-8),]

arms <- unique(sig_df$arm)

print(pheno)
for (arm in arms) {
	print(arm)
	arm_df <- merged_df[(merged_df$arm == arm),]
	out_name <- paste0(paste(pheno, arm, sep = "-"), ".txt")
	write.table(arm_df, file = paste(out_path, out_name, sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
