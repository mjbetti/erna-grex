#######################################
####Complete functions: 
####Ayush Giri
####Date: 12/09/2016

#1) new_xaxis_func: Prepares GWAS dataset to have ordered axis positions and other necessary items; 

new_xaxis_func_v2 <- function(x, chr_name="CHR", pos_name="BP", p_val_name="P", snp_name = "SNP", ...) {
  
  if ((chr_name %in% names(x)) & (pos_name %in% names(x)) & (p_val_name %in% names(x))) {
    if (is.null(x[[snp_name]])) {
      d1=data.frame(CHR=x[[chr_name]], BP=x[[pos_name]], P=x[[p_val_name]])
    } else {d1=data.frame(CHR=x[[chr_name]], BP=x[[pos_name]], P=x[[p_val_name]], SNP=x[[snp_name]])
    }
  }
  else {
    stop("Please double check columns to ensure column names for chromosome, BP and P-value match")
  }
  
  d1 <- subset(d1, is.numeric(CHR) & is.numeric(BP) & is.numeric(P))
  #sort
  d1 <- d1[order(d1$CHR, d1$BP), ]
  
  #convert to p-val to -log10 p
  d1$log_p = log10(d1$P)
  
  n_ordered_chroms <- unique(d1$CHR)
  first <- n_ordered_chroms[1]
  n_ordered_chroms[length(n_ordered_chroms)]
  
  store_last_positions = c()
  chrom_track = c()
  
  if (length(n_ordered_chroms) == 1) {
    d1$new_pos = d1$pos
  } else {
    for (chrom in n_ordered_chroms) {
      if (chrom == first) {
        d1$new_pos[d1$CHR == chrom] = d1$BP[d1$CHR == chrom]
        #last_pos <- d1$new_pos[length(d1$new_pos[d1$CHR == chrom])]
        last_pos <- as.numeric(d1$new_pos[length(d1$new_pos[!is.na(d1$new_pos)])])
        store_last_positions <- c(last_pos)
        chrom_track[chrom] <- c(paste('chromosome', chrom, sep = ' '))
      } else {
        d1$new_pos[d1$CHR == chrom] = d1$BP[d1$CHR == chrom] + last_pos
        last_pos <- as.numeric(tail(subset(d1$new_pos, d1$CHR == chrom), 1))
        store_last_positions = c(store_last_positions, last_pos)
        chrom_track[chrom] <- c(paste('chromosome', chrom, sep = ' '))
      }
      df_for_ticks <- NULL
      for (i in unique(d1$CHR)) {
        df_for_ticks$chrs[i] <- i
        df_for_ticks$chr_mid_points[which(df_for_ticks$chrs == i)] <- median(d1$new_pos[which(d1$CHR == i)])
      }
      df_for_ticks <- as.data.frame(list(df_for_ticks))      
    }
    #out_data <- data.frame(d1)
    #return(out_data)
    #result = list(d1, store_last_positions, chrom_track)
  }
  
  df_for_ticks <- cbind(df_for_ticks, store_last_positions)
  
  d1$col_cat2 <- NA
  for (i in 1:length(n_ordered_chroms)) {
    d1$col_cat[which(d1$CHR == n_ordered_chroms[i])] <- as.numeric(i)
    d1$col_cat2[which(d1$col_cat == i)] <- as.numeric(i) %% 2
  }
    #d1$col_cat2 <- ifelse(d1$col_cat %% 2 == 1, 1, 2)
  
  #result = list(d1, store_last_positions, chrom_track, n_ordered_chroms, df_for_ticks)
  result = list(d1, df_for_ticks)
  return(result)
}



#################################mirror plot modification to include bells and whistles
#2) The actual plotting function; 

     mirror_plot_func <- function(d1=d1, new_data=new_data, chr = NULL, y_min = NULL, y_max = NULL, y_ticks = NULL, zoom_left = NULL, zoom_right = NULL, gene_tag_p = NULL, x=df_for_ticks, chr_mid_points = "chr_mid_points", chrs = "chrs", store_last_positions = "store_last_positions", color_tissue = NULL, sig_line1 = NULL, sig_line1_sug = NULL, sig_line2 = NULL, sig_line2_sug = NULL, sig_line1_color = NULL, sig_line1_sug_color = NULL, sig_line2_color = NULL, sig_line2_sug_color = NULL) {
      
       df_for_ticks=data.frame(chr_mid_points=df_for_ticks[[chr_mid_points]], chrs=df_for_ticks[[chrs]], store_last_positions=df_for_ticks[[store_last_positions]])
       labels_cat <- c(sort(unique(as.character(new_data$dtype))))
       test_ylab <- expression(log["10"]*italic((p))~"&"~-log["10"]*italic((p)))
	
	#Setting default y-minima and y-maxima for plots; 
	#if(is.null(y_min) == TRUE) {y_min <- round(min(d1$log_p, na.rm=TRUE) - 1)} else {y_min <- y_min}
	#if(is.null(y_max) == TRUE) {y_max <- round(max(new_data$log_p, na.rm=TRUE) + 1)} else {y_max <- y_max}
	#if(is.null(y_ticks) == TRUE) {
		#if((y_max - y_min) > 16) {y_ticks <- 16} else {y_ticks <- round(y_max - y_min)} 
	#} else {y_ticks <- y_ticks}
	
       
	##Modifying gene-specific dataset for plotting
       new_data$new_pos <- NA
       new_data$col_cat2 <- NA
       for (i in order(unique(new_data$CHR))) {
         if(i == df_for_ticks$chrs[1]) {
           new_data$new_pos[which(new_data$CHR == i)] <- new_data$MID_POS[which(new_data$CHR == i)]
         }
         else {
           new_data$new_pos[which(new_data$CHR == i)] <- new_data$MID_POS[which(new_data$CHR == i)] + df_for_ticks$store_last_positions[which(df_for_ticks$chrs == i)-1]
         }
       }
       
       new_data$col_cat2 <- NA
       for (i in 1:length(order(unique(new_data$CHR)))) {
         new_data$col_cat[which(new_data$CHR == order(unique(new_data$CHR))[i])] <- as.numeric(i)
         new_data$col_cat2[which(new_data$col_cat == i)] <- as.numeric(i) %% 2
       }
       
  ##Color Tissue Null option      
       
       if(is.null(color_tissue) == TRUE) {
         #color_tissue =  c("firebrick", "blue", "goldenrod", "#D95F02", "#7570B3", "#E7298A", "forestgreen", "darkgray", "black", "dodgerblue")
	color_tissue =  c("bisque3","peachpuff2","lightgreen","firebrick4","red","orangered", "powderblue", "blue4","darkturquoise","dodgerblue","cyan","steelblue", 
			"lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "turquoise", "midnightblue", "lightpink1","gray32","grey85","goldenrod4","khaki3", 
			"lightgoldenrod","tan2","yellow","tomato","tomato3", "burlywood4","plum", "lawngreen", "chocolate4","paleturquoise4","palevioletred3","seagreen1",
			"darkgreen","honeydew4", "cornsilk2","wheat2","sandybrown","lightsalmon3","gold","darkslategrey","forestgreen","lightcoral", "deeppink","gray0","darkseagreen3")

       }
       else {
         color_tissue = color_tissue 
       }
       
      
       
      if((is.null(chr) & is.null(zoom_left) & is.null(zoom_right)) == TRUE) {
        
		#Setting default y-minima and y-maxima for plots; 
		if(is.null(y_min) == TRUE) {y_min <- round(min(d1$log_p, na.rm=TRUE) - 1)} else {y_min <- y_min}
		if(is.null(y_max) == TRUE) {y_max <- round(max(new_data$log_p, na.rm=TRUE) + 1)} else {y_max <- y_max}
		if(is.null(y_ticks) == TRUE) {
			if((y_max - y_min) > 16) {y_ticks <- 16} else {y_ticks <- round(y_max - y_min)} 
		} else {y_ticks <- y_ticks}
        
    break_length <- round(round(y_max - y_min)/y_ticks)
    
    for_tag <- subset(new_data, log_p > gene_tag_p)
    for_tag <- for_tag[order(for_tag$Gene, -for_tag$log_p),]
    for_tag <- for_tag[!duplicated(for_tag$Gene),]
    if(is.null(sig_line1) == TRUE) {sig_line1 <- 0} else{sig_line1 <- sig_line1}
    if(is.null(sig_line1_sug) == TRUE) {sig_line1_sug <- 0} else{sig_line1_sug <- sig_line1_sug}
    if(is.null(sig_line2) == TRUE) {sig_line2 <- 0} else{sig_line2 <- sig_line2}
    if(is.null(sig_line2_sug) == TRUE) {sig_line2_sug <- 0} else{sig_line2_sug <- sig_line2_sug}
    if(is.null(sig_line1_color) == TRUE) {sig_line1_color <- "black"} else {sig_line1_color <- sig_line1_color}
    if(is.null(sig_line1_sug_color) == TRUE) {sig_line1_sug_color <- "black"} else {sig_line1_sug_color <- sig_line1_sug_color}
    if(is.null(sig_line2_color) == TRUE) {sig_line2_color <- "black"} else {sig_line2_color <- sig_line2_color}
    if(is.null(sig_line2_sug_color) == TRUE) {sig_line2_sug_color <- "black"} else {sig_line2_sug_color <- sig_line2_sug_color}

		p1 <- ggplot() +  theme_bw() +  scale_x_continuous(breaks=c(df_for_ticks$chr_mid_points), 
                                                        labels=c(df_for_ticks$chrs)) +
       theme(axis.text.x = element_text(size=12, color = 'black'), 
             axis.text.y = element_text(size=12, color = 'black'), 
             axis.title.x = element_text(size = 12, face = "bold", color ="black"),
             axis.title.y = element_text(size = 12, face = "bold", color ="black"),
             axis.ticks.x=element_line()) +
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          #theme(legend.position="none") +
      
        geom_point(data = d1, aes(x = new_pos, y = log_p, size = 1.0, color = factor(col_cat2))) +
        xlab("Chromosomes") + ylab(test_ylab) +
       scale_y_continuous(breaks=seq(y_min, y_max, break_length)) + expand_limits(y=c(y_min, y_max)) +

       geom_point(data = new_data, aes(x = new_pos, y = log_p, size = 1.0, color = dtype, shape = factor(col_cat2))) +
        scale_colour_manual(name = "Tissue Type", values = c("darkgray", "black", color_tissue), labels = c("Single SNP", "Single SNP", labels_cat)) +
        #scale_colour_manual(values = c("darkgray", "black", color_tissue)) +
          #guides(shape = "none", size = "none", col = guide_legend(reverse = TRUE), colour = guide_legend(reverse = TRUE, override.aes = list(size=10))) + 
          guides(shape = "none", size = "none", colour = guide_legend(reverse = TRUE, override.aes = list(size=6))) + 
          
       #geom_label_repel(data = subset(new_data, log_p > gene_tag_p), aes(x = new_pos, y = log_p, label = Gene)) +
       geom_label_repel(data = for_tag, aes(x = new_pos, y = log_p, label = Gene)) +
	geom_hline(aes(yintercept = 0), size = 1) +
        geom_hline(yintercept = sig_line1, color = c(sig_line1_color), size = 1) +
        geom_hline(yintercept = sig_line1_sug, color = c(sig_line1_sug_color), size = 1) +
        geom_hline(yintercept = sig_line2, color = c(sig_line2_color), size = 1) +
        geom_hline(yintercept = sig_line2_sug, color = c(sig_line2_sug_color), size = 1)

      }
      else if((chr %in% seq(1, 25, by= 1)) & is.null(zoom_left) & is.null(zoom_right)) {
        
	d1$BP <- d1$BP/1000000
        new_data$MID_POS <- new_data$MID_POS/1000000	
		
		#Setting default y-minima and y-maxima for plots; 
		if(is.null(y_min) == TRUE) {y_min <- round(min(d1$log_p[which(d1$CHR == chr)], na.rm=TRUE) - 1)} else {y_min <- y_min}
		if(is.null(y_max) == TRUE) {y_max <- round(max(new_data$log_p[which(new_data$CHR == chr)], na.rm=TRUE) + 1)} else {y_max <- y_max}
		if(is.null(y_ticks) == TRUE) {
			if((y_max - y_min) > 16) {y_ticks <- 16} else {y_ticks <- round(y_max - y_min)} 
		} else {y_ticks <- y_ticks}
        
    break_length <- round(round(y_max - y_min)/y_ticks)

    for_tag <- subset(new_data, CHR == chr & log_p > gene_tag_p)    
    for_tag <- for_tag[order(for_tag$Gene, -for_tag$log_p),]
    for_tag <- for_tag[!duplicated(for_tag$Gene),]
   if(is.null(sig_line1) == TRUE) {sig_line1 <- 0} else{sig_line1 <- sig_line1}
   if(is.null(sig_line1_sug) == TRUE) {sig_line1_sug <- 0} else{sig_line1_sug <- sig_line1_sug}
   if(is.null(sig_line2) == TRUE) {sig_line2 <- 0} else{sig_line2 <- sig_line2}
   if(is.null(sig_line2_sug) == TRUE) {sig_line2_sug <- 0} else{sig_line2_sug <- sig_line2_sug}
    if(is.null(sig_line1_color) == TRUE) {sig_line1_color <- "black"} else {sig_line1_color <- sig_line1_color}
    if(is.null(sig_line1_sug_color) == TRUE) {sig_line1_sug_color <- "black"} else {sig_line1_sug_color <- sig_line1_sug_color}
    if(is.null(sig_line2_color) == TRUE) {sig_line2_color <- "black"} else {sig_line2_color <- sig_line2_color}
    if(is.null(sig_line2_sug_color) == TRUE) {sig_line2_sug_color <- "black"} else {sig_line2_sug_color <- sig_line2_sug_color}


        p1 <- ggplot() +  theme_bw() +  
          theme(axis.text.x = element_text(size=12, color = 'black'), 
                axis.text.y = element_text(size=12, color = 'black'), 
                axis.title.x = element_text(size = 12, face = "bold", color ="black"),
                axis.title.y = element_text(size = 12, face = "bold", color ="black"),
                axis.ticks.x=element_line()) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          #theme(legend.position="none") +
          
          geom_point(data = subset(d1, CHR == chr), aes(x = BP, y = log_p, size = 1.0, color = factor(col_cat2))) +
          xlab(paste("Chromosome", chr, "(Mb)", sep = " ")) + ylab(test_ylab) + 
	  scale_y_continuous(breaks=seq(y_min, y_max, break_length)) + expand_limits(y=c(y_min, y_max)) +
          
          geom_point(data = subset(new_data, CHR == chr), aes(x = MID_POS, y = log_p, size = 1.0, color = dtype, shape = factor(col_cat2))) +
          #scale_colour_manual(values = c("darkgray", "black", color_tissue)) +
          scale_colour_manual(name = "Tissue Type", values = c("black", color_tissue), labels = c("Single SNP", labels_cat)) +
          #scale_colour_manual(values = c("darkgray", "black", color_tissue)) +
          #guides(shape = "none", size = "none", col = guide_legend(reverse = TRUE), colour = guide_legend(reverse = TRUE, override.aes = list(size=10))) + 
          guides(shape = "none", size = "none", colour = guide_legend(reverse = TRUE, override.aes = list(size=6))) + 
         
	  #geom_label_repel(data = subset(new_data, CHR == chr & log_p > gene_tag_p), aes(x = MID_POS, y = log_p, label = Gene)) +
          geom_label_repel(data = for_tag, aes(x = MID_POS, y = log_p, label = Gene)) +

	  geom_hline(aes(yintercept = 0), size = 1) +
	#geom_hline(aes(yintercept = sig_line1), color = "red", size = 1) +
        #geom_hline(aes(yintercept = gene_tag_p), color = "red", size = 1)
        geom_hline(yintercept = sig_line1, color = c(sig_line1_color), size = 1) +
        geom_hline(yintercept = sig_line1_sug, color = c(sig_line1_sug_color), size = 1) +
        geom_hline(yintercept = sig_line2, color = c(sig_line2_color), size = 1) +
        geom_hline(yintercept = sig_line2_sug, color = c(sig_line2_sug_color), size = 1)

      }
       else if(((chr %in% seq(1, 25, by= 1)) & is.numeric(zoom_left) & is.numeric(zoom_right)) == TRUE) {
         
         ##Modifying input dataset for line graph plotting; 
         new_data1 <- data.frame(CHR=new_data$CHR, Gene=new_data$Gene, LINE_POS=new_data$START_POS, log_p=new_data$log_p, dtype=new_data$dtype, tag_identifier=0)
         new_data2 <- data.frame(CHR=new_data$CHR, Gene=new_data$Gene, LINE_POS=new_data$END_POS, log_p=new_data$log_p, dtype=new_data$dtype, tag_identifier=1)
         new_data3 <- rbind(new_data1, new_data2)
         new_data3$Gene_dtype <- paste(new_data3$Gene, new_data3$dtype, sep="_")
	 #new_data3$ENSG_dtype <- paste(new_data3$ENSG, new_data3$dtype, sep="_")
	 new_data3$LINE_POS <- new_data3$LINE_POS/1000000
         
         d1$BP <- d1$BP/1000000         

		#Setting default y-minima and y-maxima for plots; 
		if(is.null(y_min) == TRUE) {y_min <- round(min(d1$log_p[which(d1$CHR == chr & d1$BP >= zoom_left & d1$BP <= zoom_right)], na.rm=TRUE) - 1)} else {y_min <- y_min}
		if(is.null(y_max) == TRUE) {y_max <- round(max(new_data3$log_p[which(new_data3$CHR == chr & new_data3$LINE_POS >= zoom_left & new_data3$LINE_POS <= zoom_right)], na.rm=TRUE) + 1)} else {y_max <- y_max}
		if(is.null(y_ticks) == TRUE) {
			if((y_max - y_min) > 16) {y_ticks <- 16} else {y_ticks <- round(y_max - y_min)} 
			} else {y_ticks <- y_ticks}		
		 
    break_length <- round(round(y_max - y_min)/y_ticks)     
         
      for_tag <- subset(new_data3, CHR == chr & log_p > gene_tag_p & LINE_POS >= zoom_left & LINE_POS <= zoom_right) 
     for_tag <- for_tag[order(for_tag$Gene, -for_tag$log_p),]
    for_tag <- for_tag[!duplicated(for_tag$Gene),]
    if(is.null(sig_line1) == TRUE) {sig_line1 <- 0} else{sig_line1 <- sig_line1}
    if(is.null(sig_line1_sug) == TRUE) {sig_line1_sug <- 0} else{sig_line1_sug <- sig_line1_sug}
    if(is.null(sig_line2) == TRUE) {sig_line2 <- 0} else{sig_line2 <- sig_line2}
    if(is.null(sig_line2_sug) == TRUE) {sig_line2_sug <- 0} else{sig_line2_sug <- sig_line2_sug}
    if(is.null(sig_line1_color) == TRUE) {sig_line1_color <- "black"} else {sig_line1_color <- sig_line1_color}
    if(is.null(sig_line1_sug_color) == TRUE) {sig_line1_sug_color <- "black"} else {sig_line1_sug_color <- sig_line1_sug_color}
    if(is.null(sig_line2_color) == TRUE) {sig_line2_color <- "black"} else {sig_line2_color <- sig_line2_color}
    if(is.null(sig_line2_sug_color) == TRUE) {sig_line2_sug_color <- "black"} else {sig_line2_sug_color <- sig_line2_sug_color}
         
         p1 <- ggplot() +  theme_bw() +  
           theme(axis.text.x = element_text(size=12, color = 'black'), 
                 axis.text.y = element_text(size=12, color = 'black'), 
                 axis.title.x = element_text(size = 12, face = "bold", color ="black"),
                 axis.title.y = element_text(size = 12, face = "bold", color ="black"),
                 axis.ticks.x=element_line()) +
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
           #theme(legend.position="none") +
           
           geom_point(data = subset(d1, CHR == chr & BP >= zoom_left & BP <= zoom_right), aes(x = BP, y = log_p, color = factor(col_cat2))) +
           xlab(paste("Chromosome", chr, "(Mb)", sep = " ")) + ylab(test_ylab) +
           scale_y_continuous(breaks=seq(y_min, y_max, break_length)) + expand_limits(y=c(y_min, y_max)) +
           scale_x_continuous(breaks =seq(zoom_left, zoom_right, 1)) + 
           #geom_point(data = subset(new_data, CHR == chr & BP >= zoom_left & BP <= zoom_right), aes(x = BP, y = log_p, size = 1.0, color = dtype, shape = factor(col_cat2))) +
           #scale_colour_manual(values = c("darkgray", "black", color_tissue)) +
           #geom_label_repel(data = subset(new_data, CHR == chr & BP >= zoom_left & BP <= zoom_right & log_p > 1.5), aes(x = BP, y = log_p, label = Gene)) +
           #geom_hline(aes(yintercept = 0), size = 1) 
           
           #geom_line(data = subset(new_data3, CHR == chr & LINE_POS >= zoom_left & LINE_POS <= zoom_right), aes(x=LINE_POS, y = log_p, group = Gene, size = 1.0, color = dtype, shape = factor(col_cat2))) +
           geom_line(data = subset(new_data3, CHR == chr & LINE_POS >= zoom_left & LINE_POS <= zoom_right), aes(x=LINE_POS, y = log_p, group = Gene_dtype, size = 1.5, color = dtype)) +
	   #geom_line(data = subset(new_data3, CHR == chr & LINE_POS >= zoom_left & LINE_POS <= zoom_right), aes(x=LINE_POS, y = log_p, group = ENSG_dtype, size = 1.5, color = dtype)) +
           scale_colour_manual(name = "Tissue Type", values = c("black", color_tissue), labels = c("Single SNP", labels_cat)) +
           guides(shape = "none", size = "none", colour = guide_legend(reverse = TRUE, override.aes = list(size=6))) + 
           #geom_label_repel(data = subset(new_data3, CHR == chr & LINE_POS >= zoom_left & LINE_POS <= zoom_right & log_p > gene_tag_p & tag_identifier == 1), aes(x=LINE_POS, y = log_p, label = Gene)) +
           geom_label_repel(data = for_tag, aes(x = LINE_POS, y = log_p, label = Gene)) +

	   geom_hline(aes(yintercept = 0), size = 1) +
          # geom_hline(aes(yintercept = sig_line1), color = "red", size = 1) +
        #geom_hline(aes(yintercept = gene_tag_p), color = "red", size = 1)
        geom_hline(yintercept = sig_line1, color = c(sig_line1_color), size = 1) +
        geom_hline(yintercept = sig_line1_sug, color = c(sig_line1_sug_color), size = 1) +
        geom_hline(yintercept = sig_line2, color = c(sig_line2_color), size = 1) +
        geom_hline(yintercept = sig_line2_sug, color = c(sig_line2_sug_color), size = 1)

         
       }    
      else {print("If you leave CHR as blank do not enter, zoom location")}
       return(p1)
     }


