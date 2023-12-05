library("data.table")

#Declare the path of the LD scores file and TWAS effect sizes
ld_path <- "/home/bettimj/gamazon_rotation/gtex8_geno/liftover_hg19/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.pruned_geno0_maf0.05_hwe1e-6.score.ld"
twas_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all_genes/brain_only/arm/paired_erna_gene_association_loci_arm_scz_erna_and_gene_twas.txt"

#Open the LD scores as a data frame
ld_file <- fread(ld_path, header = TRUE, sep = " ", quote = "")
ld_df <- as.data.frame(ld_file)

twas_file <- fread(twas_path, header = TRUE, sep = " ", quote = "")
twas_df <- as.data.frame(twas_file)

#Brain Cerebellum
##Declare the paths of the eRNA and gene SNP weights
erna_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training/Brain_Cerebellum/weights/ENSR00000116923.txt"

