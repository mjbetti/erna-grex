source(file="/home/bettimj/aldrich_rotation/metaxcan/metaxcan_plot_functions_gtexv8.R")
library(data.table)

#############Loading GWAS file first: 

gwas_summary <- fread("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv", header= TRUE, sep = "auto")  #########Change name of GWAS file here; If needed change seperator as needed, currently it is set to whitespace(space and/or tab)

gwas_df <- new_xaxis_func_v2(gwas_summary, chr_name = "CHROM", pos_name = "POS", p_val_name = "PVAL", snp_name = "SNP") ##########Specify column names here
gwas_df_d1 <- as.data.frame(list(gwas_df[1])) ##Extracts formatted GWAS input file
df_for_ticks <- as.data.frame(list(gwas_df[2])) ##File generated from new_xaxis_func function; needed for plotting
dim(df_for_ticks)

map_file_gene <- fread("/home/bettimj/aldrich_rotation/metaxcan/gencode.v32.GRCh37.metaxcan_map.txt", sep = "auto", header = TRUE)
map_file_erna <- fread("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/ensembl87_fantom5_roadmap_hg19_reg_annos_hg19_map.txt", header = TRUE)

map_file <- rbind(map_file_gene, map_file_erna)


metaxcan_results_files <- list.files(pattern = "*.csv$")  ####################May want to change here, or move desired files to a new directory
##Create a vector with elements in the same order as your *.csv files 
c1 <- c("Brain_Amy", "Brain_Ant", "Brain_Cau", "Brain_Cerebellum", "Brain_Cer", "Brain_Cortex", "Brain_Fro", "Brain_Hippo", "Brain_Hypo", "Brain_Nuc", "Brain_Put", "Brain_Spi", "Brain_Sub", "Brain_Amy", "Brain_Ant", "Brain_Cau", "Brain_Cerebellum", "Brain_Cer", "Brain_Cortex", "Brain_Fro", "Brain_Hippo", "Brain_Hypo", "Brain_Nuc", "Brain_Put", "Brain_Spi", "Brain_Sub") #############Check to make sure the order matches metaxcan_results_files; or do ls *.csv to check order #############Check to make sure the order matches metaxcan_results_files; or do ls *.csv to check order


track_names <- data.frame(name_list = c1, file_name = unlist(metaxcan_results_files))

head(track_names)
dim(track_names)

######################Prepping MetaXcan results input files for plotting software; 

rm(appended_gene_d, a, gene_dat_generic)

###Need to complete this stuff here; 
#1) need to add the pmatch thing after the exists statement as well;
#2) Ideally you want to keep the ENSG column and use that for unique identification, then use genename column for labeling only

for (i in 1:length(track_names$name_list)) {
  
  if(!exists("appended_gene_d")){
    a <- fread(as.character(track_names$file_name[i]), header = T, sep="auto")
    
###	if(track_names$file_name[i] == "DGN-WB-unscaled0_0.5.csv") {
###	appended_gene_d <- merge(a, map_file, by.x = "gene",  by.y ="Gene")
###    names(appended_gene_d)[names(appended_gene_d) == "gene"] <- "Gene"
###	t1 <- which(duplicated(appended_gene_d$Gene))
###
###
####	if(length(t1) == 0) {appended_gene_d <- appended_gene_d} else {appended_gene_d <- appended_gene_d[-c(t1), ]}
###	} else {                                                         
    appended_gene_d <- merge(a, map_file, by.x = "gene", by.y = "ENSG")
    names(appended_gene_d)[names(appended_gene_d) == "gene"] <- "ENSG"
	t1 <- which(duplicated(appended_gene_d$ENSG))
	if(length(t1) == 0) {appended_gene_d <- appended_gene_d} else {appended_gene_d <- appended_gene_d[-c(t1), ]}
###	}
    appended_gene_d$MID_POS <- (appended_gene_d$END_POS + appended_gene_d$START_POS)/2
    appended_gene_d$log_p <- -log10(appended_gene_d$pvalue) ########## highly optional: if metaxcan p-value column name changes, then update $pvalue as needed  
    appended_gene_d$dtype <- track_names$name_list[i]
    appended_gene_d$tiss_number <- i
    rm(a)
} else if(exists("appended_gene_d")){
    a <- read.csv(file=as.character(track_names$file_name[i]), header = T)
    
###     if(track_names$file_name[i] == "DGN-WB-unscaled0_0.5.csv") {
###     gene_dat_generic <- merge(a, map_file, by.x = "gene",  by.y ="Gene")
###     names(gene_dat_generic)[names(gene_dat_generic) == "gene"] <- "Gene"
###	t1 <- which(duplicated(gene_dat_generic$Gene))
###	if(length(t1) == 0) {gene_dat_generic <- gene_dat_generic} else {gene_dat_generic <- gene_dat_generic[-c(t1), ]}
###	} else {                                                         
    gene_dat_generic <- merge(a, map_file, by.x = "gene", by.y = "ENSG")
    names(gene_dat_generic)[names(gene_dat_generic) == "gene"] <- "ENSG"
	t1 <- which(duplicated(gene_dat_generic$ENSG))
	if(length(t1) == 0) {gene_dat_generic <- gene_dat_generic} else {gene_dat_generic <- gene_dat_generic[-c(t1), ]}
###	}
    gene_dat_generic$MID_POS <- (gene_dat_generic$END_POS + gene_dat_generic$START_POS)/2
    gene_dat_generic$log_p <- -log10(gene_dat_generic$pvalue) ########## highly optional: if metaxcan p-value column name changes, then update $pvalue as needed 
    gene_dat_generic$dtype <- track_names$name_list[i]
    gene_dat_generic$tiss_number <- i
    appended_gene_d <- rbind(appended_gene_d, gene_dat_generic)
    rm(gene_dat_generic, a)
  }
  
}

appended_gene_d <- appended_gene_d[with(appended_gene_d, order(dtype, CHR, START_POS)),]
head(appended_gene_d)

########################PLOTTING BEGINS BELOW#############################
##########################################################################

library(ggplot2)
library(ggrepel)
options(bitmapType='cairo')
##Whole chromosome plot
###Optional flags gene_tag_p --- log10(p-value) threshold for naming genes
###Option color (if you have two tissues and want to specify custom colors; do: color_tissue = c("colorname1", "colorname2") and so forth. 

# tiff(filename = "Entire_Genome_result_schizophrenia_erna_jti_brain_tissues.tiff", height= 10, width = 20, unit = "in", res=300, pointsize=11, compression ='lzw') ###Optional, filename change for output if desired. 
# mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = NULL, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 100, 
# x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
# 			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
# y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 5.113, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
# dev.off()


tiff(filename = "Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", height= 10, width = 20, unit = "in", res=300, pointsize=11, compression ='lzw') ###Optional, filename change for output if desired. 
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = NULL, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

##Chromosome specific plots Note that chr = some number not NULL
###Optional flags gene_tag_p --- log10(p-value) threshold for naming genes
###Option color (if you have two tissues and want to specify custom colors; do: color_tissue = c("colorname1", "colorname2") and so forth. 

#tiff(filename = "Chr_6_v10302018_v2.tiff", height= 10, width = 20, unit = "in", res=300, pointsize=11, compression ='lzw')
#mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 6, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 6.6, 
#x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = NULL, y_min = NULL, 
#y_max = NULL, y_ticks = NULL, sig_line1 = 6.6, sig_line2 = -7.4, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
#dev.off()

###############Plot all chromosome specific plots using loop. 
# plot_list <- list()
# for (i in 1:length(unique(df_for_ticks$chrs))) {
# fname_print <- paste("Chr_", df_for_ticks$chrs[i], "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
# tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
# mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = df_for_ticks$chrs[i], zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
# x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
# 			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
# y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
# dev.off()
# }
# rm(plot_list)

fname_print <- paste("Chr_", 1, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 1, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 2, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 2, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 3, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 3, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 4, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 4, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 5, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 5, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 6, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 6, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 7, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 7, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 8, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 8, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 9, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 9, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 10, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 10, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 11, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 11, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 12, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 12, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 13, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 13, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 14, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 14, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 15, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 15, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 16, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 16, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 17, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 17, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 18, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 18, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 19, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 19, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 20, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 20, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 21, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 21, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

fname_print <- paste("Chr_", 22, "_Entire_Genome_labeled_result_schizophrenia_gene_brain_tissues_att4.tiff", sep = "") ###Modify name of file as you wish
tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 22, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.198, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
dev.off()

###Zoomed in regional plots; Note that chr = some number NOT NULL, and zoom_left and zoom_right also have numbers in Mega bases. 
###Optional flags gene_tag_p --- log10(p-value) threshold for naming genes 
###Option color (if you have two tissues and want to specify custom colors; do: color_tissue = c("colorname1", "colorname2") and so forth. 

#tiff(filename = "CASZ1.tiff", height= 10, width = 18, unit = "in", res=300, pointsize=11, compression ='lzw')
#mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 1, zoom_left = 9.7 ,  zoom_right = 11.8, 
#gene_tag_p = 6.6, x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = NULL, 
#y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.6, sig_line2 = -7.4, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
#dev.off()

#tiff(filename = "WNT2B.tiff", height= 10, width = 18, unit = "in", res=300, pointsize=11, compression ='lzw')
#mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 1, zoom_left = 112 ,  zoom_right = 114, 
#gene_tag_p = 6.6, x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = NULL, 
#y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.6, sig_line2 = -7.4, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
#dev.off()

#tiff(filename = "KCNK3.tiff", height= 10, width = 18, unit = "in", res=300, pointsize=11, compression ='lzw')
#mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 2, zoom_left = 25.9 ,  zoom_right = 27.9, 
#gene_tag_p = 6.6, x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = NULL, 
#y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.6, sig_line2 = -7.4, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
#dev.off()

#tiff(filename = "FBN2.tiff", height= 10, width = 18, unit = "in", res=300, pointsize=11, compression ='lzw')
#mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 5, zoom_left = 126.6 ,  zoom_right = 128.8, 
#gene_tag_p = 6.6, x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = NULL, 
#y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.6, sig_line2 = -7.4, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
#dev.off()

#tiff(filename = "HOXA13.tiff", height= 10, width = 18, unit = "in", res=300, pointsize=11, compression ='lzw')
#mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 7, zoom_left = 26.2 ,  zoom_right = 28.2, 
#gene_tag_p = 6.6, x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = NULL, 
#y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.6, sig_line2 = -7.4, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
#dev.off()

#tiff(filename = "BLK.tiff", height= 10, width = 18, unit = "in", res=300, pointsize=11, compression ='lzw')
#mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 8, zoom_left = 10.3,  zoom_right = 12.4, 
#gene_tag_p = 6.6, x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = NULL, 
#y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.6, sig_line2 = -7.4, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
#dev.off()

#tiff(filename = "LSP1.tiff", height= 10, width = 18, unit = "in", res=300, pointsize=11, compression ='lzw')
#mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 11, zoom_left = 0.8,  zoom_right = 2.9, 
#gene_tag_p = 6.6, x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = NULL, 
#y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.6, sig_line2 = -7.4, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
#dev.off()

#tiff(filename = "RXFP2.tiff", height= 10, width = 18, unit = "in", res=300, pointsize=11, compression ='lzw')
#mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 13, zoom_left = 31.3,  zoom_right = 33.4, 
#gene_tag_p = 6.6, x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = NULL, 
#y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.6, sig_line2 = -7.4, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
#dev.off()

#tiff(filename = "Chr13_intergenic.tiff", height= 10, width = 18, unit = "in", res=300, pointsize=11, compression ='lzw')
#mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 13, zoom_left = 72.8 ,  zoom_right = 73.8, 
#gene_tag_p = 6.6, x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = NULL, 
#y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.6, sig_line2 = -7.4, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
#dev.off()

#tiff(filename = "Chr12_intergenic.tiff", height= 10, width = 18, unit = "in", res=300, pointsize=11, compression ='lzw')
#mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 12, zoom_left = 114.5 ,  zoom_right = 116.5, 
#gene_tag_p = 6.6, x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = NULL, 
#y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.6, sig_line2 = -7.4, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
#dev.off()

#tiff(filename = "Chr6_intergenic.tiff", height= 10, width = 18, unit = "in", res=300, pointsize=11, compression ='lzw')
#mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 6, zoom_left = 126 ,  zoom_right = 128, 
#gene_tag_p = 6.6, x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = NULL, 
#y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.6, sig_line2 = -7.4, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
#dev.off()

