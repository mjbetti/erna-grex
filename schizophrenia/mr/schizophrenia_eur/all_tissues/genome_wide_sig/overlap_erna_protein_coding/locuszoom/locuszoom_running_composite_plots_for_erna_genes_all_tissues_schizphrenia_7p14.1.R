source(file="/home/bettimj/aldrich_rotation/metaxcan/miami_metaxcan_plot_functions_gtexv8_yaxis_fixed.R")
library(data.table)

#############Loading GWAS file first: 

gwas_summary <- fread("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv", header= TRUE, sep = "auto")  #########Change name of GWAS file here; If needed change seperator as needed, currently it is set to whitespace(space and/or tab)

gwas_df <- new_xaxis_func_v2(gwas_summary, chr_name = "CHROM", pos_name = "POS", p_val_name = "PVAL", snp_name = "SNP") ##########Specify column names here
gwas_df_d1 <- as.data.frame(list(gwas_df[1])) ##Extracts formatted GWAS input file
df_for_ticks <- as.data.frame(list(gwas_df[2])) ##File generated from new_xaxis_func function; needed for plotting
dim(df_for_ticks)

map_file_gene <- fread("/home/bettimj/aldrich_rotation/metaxcan/gencode.v32.GRCh37.metaxcan_map.txt", sep = "auto", header = TRUE)
map_file_erna <- fread("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/ensembl87_fantom5_roadmap_hg19_reg_annos_hg19_map.txt", header = TRUE)

metaxcan_results_files_gene <- list.files(path = "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all_genes", pattern = "schizophrenia2022_*")  ####################May want to change here, or move desired files to a new directory
metaxcan_results_files_erna <- list.files(path = "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all", pattern = "schizophrenia2022_*")

##Create a vector with elements in the same order as your *.csv files 
c1 <- c("Adipose Subcutaneous", "Adipose Visceral Omentum", "Adrenal Gland", "Artery Aorta", "Artery Coronary", "Artery Tibial", "Brain Amygdala", "Brain Anterior cingulate cortex BA24", "Brain Caudate basal ganglia", "Brain Cerebellar Hemisphere", "Brain Cerebellum", "Brain Cortex", "Brain Frontal Cortex_BA9", "Brain Hippocampus", "Brain Hypothalamus", "Brain Nucleus accumbens basal ganglia", "Brain Putamen basal ganglia", "Brain Spinal cord cervical c-1", "Brain Substantia nigra", "Breast Mammary Tissue", "Cells Cultured fibroblasts", "Cells EBV-transformed lymphocytes", "Colon Sigmoid", "Colon Transverse", "Esophagus Gastroesophageal Junction", "Esophagus Mucosa", "Esophagus Muscularis", "Heart Atrial Appendage", "Heart Left Ventricle", "Kidney Cortex", "Liver", "Lung", "Minor Salivary Gland", "Muscle Skeletal", "Nerve Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin Not Sun Exposed Suprapubic", "Skin Sun Exposed Lower leg", "Small Intestine Terminal Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole Blood") #############Check to make sure the order matches metaxcan_results_files; or do ls *.csv to check order #############Check to make sure the order matches metaxcan_results_files; or do ls *.csv to check order

track_names_gene <- data.frame(name_list = c1, file_name = unlist(metaxcan_results_files_gene))
track_names_erna <- data.frame(name_list = c1, file_name = unlist(metaxcan_results_files_erna))

head(track_names_gene)
dim(track_names_gene)
head(track_names_erna)
dim(track_names_erna)

######################Prepping MetaXcan results input files for plotting software; 

rm(appended_gene_d_gene, a_gene, gene_dat_generic_gene)
rm(appended_gene_d_erna, a_erna, gene_dat_generic_erna)

###Need to complete this stuff here; 
#1) need to add the pmatch thing after the exists statement as well;
#2) Ideally you want to keep the ENSG column and use that for unique identification, then use genename column for labeling only

#####Gene
for (i in 1:length(track_names_gene$name_list)) {
  
  if(!exists("appended_gene_d_gene")){
    a_gene <- fread(paste0("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all_genes/", as.character(track_names_gene$file_name[i])), header = T, sep="auto")
    
###	if(track_names$file_name[i] == "DGN-WB-unscaled0_0.5.csv") {
###	appended_gene_d <- merge(a, map_file, by.x = "gene",  by.y ="Gene")
###    names(appended_gene_d)[names(appended_gene_d) == "gene"] <- "Gene"
###	t1 <- which(duplicated(appended_gene_d$Gene))
###
###
####	if(length(t1) == 0) {appended_gene_d <- appended_gene_d} else {appended_gene_d <- appended_gene_d[-c(t1), ]}
###	} else {                                                         
    appended_gene_d_gene <- merge(a_gene, map_file_gene, by.x = "gene", by.y = "ENSG")
    names(appended_gene_d_gene)[names(appended_gene_d_gene) == "gene"] <- "ENSG"
	t1_gene <- which(duplicated(appended_gene_d_gene$ENSG))
	if(length(t1_gene) == 0) {appended_gene_d_gene <- appended_gene_d_gene} else {appended_gene_d_gene <- appended_gene_d_gene[-c(t1_gene), ]}
###	}
    appended_gene_d_gene$MID_POS <- (appended_gene_d_gene$END_POS + appended_gene_d_gene$START_POS)/2
    appended_gene_d_gene$log_p <- -log10(appended_gene_d_gene$pvalue) ########## highly optional: if metaxcan p-value column name changes, then update $pvalue as needed  
    appended_gene_d_gene$dtype <- track_names_gene$name_list[i]
    appended_gene_d_gene$tiss_number <- i
    rm(a_gene)
} else if(exists("appended_gene_d_gene")){
    a_gene <- read.csv(file = paste0("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all_genes/", as.character(track_names_gene$file_name[i])), header = T)
    
###     if(track_names$file_name[i] == "DGN-WB-unscaled0_0.5.csv") {
###     gene_dat_generic <- merge(a, map_file, by.x = "gene",  by.y ="Gene")
###     names(gene_dat_generic)[names(gene_dat_generic) == "gene"] <- "Gene"
###	t1 <- which(duplicated(gene_dat_generic$Gene))
###	if(length(t1) == 0) {gene_dat_generic <- gene_dat_generic} else {gene_dat_generic <- gene_dat_generic[-c(t1), ]}
###	} else {                                                         
    gene_dat_generic_gene <- merge(a_gene, map_file_gene, by.x = "gene", by.y = "ENSG")
    names(gene_dat_generic_gene)[names(gene_dat_generic_gene) == "gene"] <- "ENSG"
	t1_gene <- which(duplicated(gene_dat_generic_gene$ENSG))
	if(length(t1_gene) == 0) {gene_dat_generic_gene <- gene_dat_generic_gene} else {gene_dat_generic_gene <- gene_dat_generic_gene[-c(t1_gene), ]}
###	}
    gene_dat_generic_gene$MID_POS <- (gene_dat_generic_gene$END_POS + gene_dat_generic_gene$START_POS)/2
    gene_dat_generic_gene$log_p <- -log10(gene_dat_generic_gene$pvalue) ########## highly optional: if metaxcan p-value column name changes, then update $pvalue as needed 
    gene_dat_generic_gene$dtype <- track_names_gene$name_list[i]
    gene_dat_generic_gene$tiss_number <- i
    appended_gene_d_gene <- rbind(appended_gene_d_gene, gene_dat_generic_gene)
    rm(gene_dat_generic_gene, a_gene)
  }
  
}

appended_gene_d_gene <- appended_gene_d_gene[with(appended_gene_d_gene, order(dtype, CHR, START_POS)),]
head(appended_gene_d_gene)

#####eRNA
for (i in 1:length(track_names_erna$name_list)) {
  
  if(!exists("appended_gene_d_erna")){
    a_erna <- fread(paste0("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all/", as.character(track_names_erna$file_name[i])), header = T, sep="auto")
    
###	if(track_names$file_name[i] == "DGN-WB-unscaled0_0.5.csv") {
###	appended_gene_d <- merge(a, map_file, by.x = "gene",  by.y ="Gene")
###    names(appended_gene_d)[names(appended_gene_d) == "gene"] <- "Gene"
###	t1 <- which(duplicated(appended_gene_d$Gene))
###
###
####	if(length(t1) == 0) {appended_gene_d <- appended_gene_d} else {appended_gene_d <- appended_gene_d[-c(t1), ]}
###	} else {                                                         
    appended_gene_d_erna <- merge(a_erna, map_file_erna, by.x = "gene", by.y = "ENSG")
    names(appended_gene_d_erna)[names(appended_gene_d_erna) == "gene"] <- "ENSG"
	t1_erna <- which(duplicated(appended_gene_d_erna$ENSG))
	if(length(t1_erna) == 0) {appended_gene_d_erna <- appended_gene_d_erna} else {appended_gene_d_erna <- appended_gene_d_erna[-c(t1_erna), ]}
###	}
    appended_gene_d_erna$MID_POS <- (appended_gene_d_erna$END_POS + appended_gene_d_erna$START_POS)/2
    appended_gene_d_erna$log_p <- -log10(appended_gene_d_erna$pvalue) ########## highly optional: if metaxcan p-value column name changes, then update $pvalue as needed  
    appended_gene_d_erna$dtype <- track_names_erna$name_list[i]
    appended_gene_d_erna$tiss_number <- i
    rm(a_erna)
} else if(exists("appended_gene_d_erna")){
    a_erna <- read.csv(file = paste0("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all/", as.character(track_names_erna$file_name[i])), header = T)
    
###     if(track_names$file_name[i] == "DGN-WB-unscaled0_0.5.csv") {
###     gene_dat_generic <- merge(a, map_file, by.x = "gene",  by.y ="Gene")
###     names(gene_dat_generic)[names(gene_dat_generic) == "gene"] <- "Gene"
###	t1 <- which(duplicated(gene_dat_generic$Gene))
###	if(length(t1) == 0) {gene_dat_generic <- gene_dat_generic} else {gene_dat_generic <- gene_dat_generic[-c(t1), ]}
###	} else {                                                         
    gene_dat_generic_erna <- merge(a_erna, map_file_erna, by.x = "gene", by.y = "ENSG")
    names(gene_dat_generic_erna)[names(gene_dat_generic_erna) == "gene"] <- "ENSG"
	t1_erna <- which(duplicated(gene_dat_generic_erna$ENSG))
	if(length(t1_erna) == 0) {gene_dat_generic_erna <- gene_dat_generic_erna} else {gene_dat_generic_erna <- gene_dat_generic_erna[-c(t1_erna), ]}
###	}
    gene_dat_generic_erna$MID_POS <- (gene_dat_generic_erna$END_POS + gene_dat_generic_erna$START_POS)/2
    gene_dat_generic_erna$log_p <- -log10(gene_dat_generic_erna$pvalue) ########## highly optional: if metaxcan p-value column name changes, then update $pvalue as needed 
    gene_dat_generic_erna$dtype <- track_names_erna$name_list[i]
    gene_dat_generic_erna$tiss_number <- i
    appended_gene_d_erna <- rbind(appended_gene_d_erna, gene_dat_generic_erna)
    rm(gene_dat_generic_erna, a_erna)
  }
  
}

appended_gene_d_erna <- appended_gene_d_erna[with(appended_gene_d_erna, order(dtype, CHR, START_POS)),]
head(appended_gene_d_erna)

########################PLOTTING BEGINS BELOW#############################
##########################################################################

library(ggplot2)
library(ggrepel)
options(bitmapType='cairo')
##Whole chromosome plot
###Optional flags gene_tag_p --- log10(p-value) threshold for naming genes
###Option color (if you have two tissues and want to specify custom colors; do: color_tissue = c("colorname1", "colorname2") and so forth. 

# tiff(filename = "miami_Entire_Genome_result_lung_cancer_erna_gene_1.36e-06_2.5e-07_100_all_tissues.tiff", height= 10, width = 20, unit = "in", res=300, pointsize=11, compression ='lzw') ###Optional, filename change for output if desired. 
# mirror_plot_func(d1=appended_gene_d_gene, new_data=appended_gene_d_erna, chr = NULL, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 100, 
# x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
# 			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue"), 
# y_min = -24, y_max = 24, y_ticks = NULL, sig_line1 = 6.157, sig_line2 = -5.113, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
# dev.off()


# tiff(filename = "miami_Entire_Genome_result_schizophrenia_erna_gene_all_tissues.tiff", height= 10, width = 20, unit = "in", res=300, pointsize=11, compression ='lzw') ###Optional, filename change for output if desired. 
# mirror_plot_func(d1=appended_gene_d_erna, new_data=appended_gene_d_gene, chr = NULL, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
# x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("bisque3","peachpuff2","lightgreen","firebrick4","red","orangered", "powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
# 			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue", "lightpink1","gray32","grey85","goldenrod4","khaki3", 
# 			"lightgoldenrod","tan2","yellow","tomato","tomato3", "burlywood4","plum", "lawngreen", "chocolate4","paleturquoise4","palevioletred3","seagreen1",
# 			"darkgreen","honeydew4", "cornsilk2","wheat2","sandybrown","lightsalmon3","gold","darkslategrey","forestgreen","lightcoral", "deeppink","gray0","darkseagreen3"), 
# y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 5.911, sig_line1_sug = 5.461, sig_line2 = -6.838, sig_line2_sug = -5.718, sig_line1_color =  c("darkred"), sig_line1_sug_color =  c("blue"), sig_line2_color = c("darkred"), sig_line2_sug_color = c("blue"))
# dev.off()

##Chromosome specific plots Note that chr = some number not NULL
###Optional flags gene_tag_p --- log10(p-value) threshold for naming genes
###Option color (if you have two tissues and want to specify custom colors; do: color_tissue = c("colorname1", "colorname2") and so forth. 

#tiff(filename = "Chr_6_v10302018_v2.tiff", height= 10, width = 20, unit = "in", res=300, pointsize=11, compression ='lzw')
#mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = 6, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 6.6, 
#x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = NULL, y_min = NULL, 
#y_max = NULL, y_ticks = NULL, sig_line1 = 6.6, sig_line2 = -7.4, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
#dev.off()

###############Plot all chromosome specific plots using loop. 
#plot_list <- list()
#for (i in 1:length(unique(df_for_ticks$chrs))) {
#fname_print <- paste("Chr_", df_for_ticks$chrs[i], "_v10092017.tiff", sep = "") ###Modify name of file as you wish
#tiff(filename = fname_print, height = 10, width = 20, unit ="in", res =300, pointsize =11, compression = 'lzw')
#plot <- mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = df_for_ticks$chrs[i], zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 4, 
#x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = NULL, 
#y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.6, sig_line2 = -7.4, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))

#plot_list[[i]] = plot
#print(plot_list[[i]])
#dev.off()
#}
#rm(plot_list)


###Zoomed in regional plots; Note that chr = some number NOT NULL, and zoom_left and zoom_right also have numbers in Mega bases. 
###Optional flags gene_tag_p --- log10(p-value) threshold for naming genes 
###Option color (if you have two tissues and want to specify custom colors; do: color_tissue = c("colorname1", "colorname2") and so forth. 

tiff(filename = "7p14.1_miami_Entire_Genome_result_schizophrenia_erna_gene_all_tissues.tiff", height= 10, width = 18, unit = "in", res=300, pointsize=11, compression ='lzw')
mirror_plot_func(d1=appended_gene_d_erna, new_data=appended_gene_d_gene, chr = 7, zoom_left = 37.2,  zoom_right = 43.3, 
gene_tag_p = 3.0, x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("bisque3","peachpuff2","lightgreen","firebrick4","red","orangered", "powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
 			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue", "lightpink1","gray32","grey85","goldenrod4","khaki3", 
 			"lightgoldenrod","tan2","yellow","tomato","tomato3", "burlywood4","plum", "lawngreen", "chocolate4","paleturquoise4","palevioletred3","seagreen1",
 			"darkgreen","honeydew4", "cornsilk2","wheat2","sandybrown","lightsalmon3","gold","darkslategrey","forestgreen","lightcoral", "deeppink","gray0","darkseagreen3"), 
y_min = -22, y_max = 16, y_ticks = NULL, sig_line1 = 5.911, sig_line1_sug = 5.461, sig_line2 = -6.838, sig_line2_sug = -5.718, sig_line1_color =  c("darkred"), sig_line1_sug_color =  c("blue"), sig_line2_color = c("darkred"), sig_line2_sug_color = c("blue"))
dev.off()

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

