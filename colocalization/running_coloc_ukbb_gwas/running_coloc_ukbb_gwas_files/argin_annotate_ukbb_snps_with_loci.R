library("dplyr")
library("data.table")
library("optparse")

option_list = list(
	make_option(c("-i", "--in_path"), type = "character", default = NULL, help = "path of input summary stats file", metavar = "character"),
	make_option(c("-l", "--loci_path"), type = "character", default = NULL, help = "path of file containing chromosomal loci", metavar = "character"),
	make_option(c("-c", "--chromosome"), type = "numeric", default = NULL, help = "chromosome number", metavar = "numeric"),
	make_option(c("-o", "--output_path"), type = "character", default = NULL, help = "path of the output directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#Declare the output directory to which these results will be saved
out_dir <- opt$output_path

#Declare the path of the AALC GWAS results file
in_path <- opt$in_path

#Declare the path of the hg38 cytoBand file downloaded from UCSC. This will be used to label SNPs based on their chromosomal location (arm).
cytoband_path <- opt$loci_path

#Open the AALC concatenated output data as a data frame, and assign the p-values to a vector
in_file <- fread(in_path, header = TRUE, sep = "\t", quote = "")
in_df <- as.data.frame(in_file)
#p_vals <- in_df$P

#in_df_sig <- na.omit(in_df[(in_df$P <= 1e-6),])
#in_df_nosig <- na.omit(in_df[(in_df$P > 1e-6),])
in_df_sig <- in_df

#Open the cytoBand file as a data frame
cytoband_file <- fread(cytoband_path, header = FALSE, sep = "\t", quote = "")
cytoband_df <- as.data.frame(cytoband_file)

#Ignore
#################################################################
#Define the function that we will use to match up SNP coordinates with their cytoBand range to return chromosome arm (adapted from https://stackoverflow.com/questions/24766104/checking-if-value-in-vector-is-in-range-of-values-in-different-length-vector)
getValue <- function(x, data) {
  tmp <- data %>%
    filter(V2 <= x, x <= V3)
  return(tmp$V4)
}

#Loop through and split up the coordinates by chromosome
in_df_sig_chr = in_df_sig[(in_df_sig$CHROM == opt$chromosome),]
cytoband_df_chr = cytoband_df[(cytoband_df[,1] == paste0("chr", opt$chromosome)),]
in_snp_coord = in_df_sig_chr$POS
arm_vector = sapply(in_snp_coord, getValue, cytoband_df_chr)
arm_vector = paste0(opt$chromosome, arm_vector)
in_df_sig_chr = cbind(in_df_sig_chr, arm_vector)

#Output our full, concatenated data frame with chromosomal arm annotations to a new file
final_out_name <- paste(paste0("chr", opt$chromosome), "ukbb_snps.arm.txt", sep = "_")
write.table(in_df_sig_chr, file = paste(opt$output_path, final_out_name, sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)