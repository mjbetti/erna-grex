library("optparse")
library("data.table")

option_list = list(
  make_option("--input", action="store", default=NA, type='character',
              help="input file path [required]"),
  make_option("--in_dir", action="store", default=NA, type='character',
              help="the directory containing the input file"),
  make_option("--out_dir", action="store", default=NA, type='character',
              help="the directory to which the output will be saved")
)

opt = parse_args(OptionParser(option_list=option_list))

#Open the input file as a data frame
in_file <- fread(paste(opt$in_dir, opt$input, sep = "/"), header = TRUE, sep = "\t", quote = "")
in_df <- as.data.frame(in_file)

#Split the variant ID column up and reformat in the same format as our eRNA models
split_df <- strsplit(in_df[,1], "[:]")
chr <- unlist(lapply(split_df, `[`, 1))
pos <- unlist(lapply(split_df, `[`, 2))
maj_allele <- unlist(lapply(split_df, `[`, 3))
min_allele <- unlist(lapply(split_df, `[`, 4))

new_varid <- paste(paste0("chr", chr), pos, maj_allele, min_allele, "b37", sep = "_")
in_df$variant <- new_varid

#Write out the output
write.table(in_df, file = paste(opt$out_dir, paste0(opt$input, ".erna_reformat.tsv"), sep = "/"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)