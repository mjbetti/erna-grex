source(file="/home/bettimj/aldrich_rotation/metaxcan/metaxcan_plot_functions_gtexv8.R")
library(data.table)

#############Loading GWAS file first: 

gwas_summary <- fread("/home/bettimj/gamazon_rotation/gwas_catalog/schizophrenia/for_metaxcan_in/snpid_daner_PGC_SCZ52_0513a.hq2", header= TRUE, sep = "auto")  #########Change name of GWAS file here; If needed change seperator as needed, currently it is set to whitespace(space and/or tab)

gwas_df <- new_xaxis_func_v2(gwas_summary, chr_name = "CHR", pos_name = "BP", p_val_name = "P", snp_name = "SNP") ##########Specify column names here
gwas_df_d1 <- as.data.frame(list(gwas_df[1])) ##Extracts formatted GWAS input file
df_for_ticks <- as.data.frame(list(gwas_df[2])) ##File generated from new_xaxis_func function; needed for plotting
dim(df_for_ticks)

map_file <- fread("/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/ensembl87_fantom5_roadmap_hg19_reg_annos_hg19_map.txt", sep = "auto", header = TRUE)


metaxcan_results_files <- list.files(pattern = "*.csv$")  ####################May want to change here, or move desired files to a new directory
##Create a vector with elements in the same order as your *.csv files 
c1 <- c("Adipose_sub", "Adipose_vis", "Adrenal", "Artery_aor", "Artery_cor", "Artery_tib", "Brain_Amy", "Brain_Ant", "Brain_Cau", "Brain_Cerebellum", "Brain_Cer", "Brain_Cortex", "Brain_Fro", "Brain_Hippo", "Brain_Hypo", "Brain_Nuc", "Brain_Put", "Brain_Spi", "Brain_Sub", "Breast", "Cells_ebv", "Cells_fibroblast", "Colon_sigmoid", "Colon_transverse", "Esophagus_gast", "Esophagus_muc", "Esophagus_mus", "Heart_Atr", "Heart_Left", "Kidney", "Liver", "Lung", "Minor_salivary_gland", "Muscle_skeletal", "Nerve", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_no_sun", "Skin_sun", "Small_intestine", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_blood") #############Check to make sure the order matches metaxcan_results_files; or do ls *.csv to check order #############Check to make sure the order matches metaxcan_results_files; or do ls *.csv to check order


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


tiff(filename = "Entire_Genome_labeled_result_schizophrenia_gene_all_tissues_att4.tiff", height= 10, width = 20, unit = "in", res=300, pointsize=11, compression ='lzw') ###Optional, filename change for output if desired. 
mirror_plot_func(d1=gwas_df_d1, new_data=appended_gene_d, chr = NULL, zoom_left = NULL,  zoom_right = NULL, gene_tag_p = 3.0, 
x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", color_tissue = c("bisque3","peachpuff2","lightgreen","firebrick4","red","orangered", "salmon4", "powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue", "lightpink1","gray32","grey85","palegreen1", "palegreen2", "goldenrod4","khaki3", 
			"lightgoldenrod","tan2","yellow", "deeppink4", "tomato","tomato3", "burlywood4","plum", "lawngreen", "chocolate4","paleturquoise4","palevioletred3","seagreen1",
			"darkgreen","honeydew4", "cornsilk2","wheat2","sandybrown","lightsalmon3","gold","darkslategrey","forestgreen","lightcoral", "deeppink","gray0","darkseagreen3"), 
y_min = NULL, y_max = NULL, y_ticks = NULL, sig_line1 = 6.838, sig_line2 = -7.30, sig_line1_color =  c("darkred"), sig_line2_color = c("darkred"))
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

