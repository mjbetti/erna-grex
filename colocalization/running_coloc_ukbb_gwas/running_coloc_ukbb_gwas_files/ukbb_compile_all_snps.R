library("data.table")

in_dir <- "/data/g_gamazon_lab/UKBioBank/both_sexes"
out_dir <- "/home/bettimj/gamazon_rotation/ukbb_with_arm"

file <- fread(paste(in_dir, "20516.gwas.imputed_v3.both_sexes.tsv", sep = "/"), quote = FALSE, sep = "\t", header = TRUE)

df <- as.data.frame(file)

df <- df[,c("variant", "rsid")]

snp <- strsplit(df$variant, ":")
chr <- unlist(lapply(snp, `[`, 1))
pos <- unlist(lapply(snp, `[`, 2))
ref <- unlist(lapply(snp, `[`, 3))
alt <- unlist(lapply(snp, `[`, 4))

varid <- paste(paste0("chr", chr), pos, ref, alt, "b37", sep = "_")

df$varid <- varid
df$chr <- chr
df$pos <- pos

names(df)[4:5] <- c("CHROM", "POS")

write.table(df, file = paste(out_dir, "ukbb_snps.txt", sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)