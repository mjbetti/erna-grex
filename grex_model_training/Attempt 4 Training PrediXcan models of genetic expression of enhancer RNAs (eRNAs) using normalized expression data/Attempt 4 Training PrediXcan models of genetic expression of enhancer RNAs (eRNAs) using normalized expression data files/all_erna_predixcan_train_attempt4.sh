module load PLINK/1.9b_5.2

#Declare the path of the main output directory, GTEx v8 genotype and expression datasets, as well as the directory of the PrediXcan script
MAIN_DIR=/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training
GENO_DIR=/home/bettimj/gamazon_rotation/gtex8_geno/liftover_hg19
EXP_DIR=/home/bettimj/gamazon_rotation/herna/normalized_to_covs
PREDIXCAN_DIR=/home/bettimj/MR-JTI/model_training/predixcan/src
ANNO_FILE=/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/ensembl87_fantom5_roadmap_hg19_reg_annos_hg19.txt

TISSUES=("Adipose_sub" "Adipose_vis" "Adrenal" "Artery_aor" "Artery_cor" "Artery_tib" "Brain_Amy" "Brain_Ant" "Brain_Cau" "Brain_Cerebellum" "Brain_Cer" "Brain_Cortex" "Brain_Fro" "Brain_Hippo" "Brain_Hypo" "Brain_Nuc" "Brain_Put" "Brain_Spi" "Brain_Sub" "Breast" "Cells_ebv" "Cells_fibroblast" "Colon_sigmoid" "Colon_transverse" "Esophagus_gast" "Esophagus_muc" "Esophagus_mus" "Heart_Atr" "Heart_Left" "Kidney" "Liver" "Lung" "Minor_salivary_gland" "Muscle_skeletal" "Nerve" "Ovary" "Pancreas" "Pituitary" "Prostate" "Skin_no_sun" "Skin_sun" "Small_intestine" "Spleen" "Stomach" "Testis" "Thyroid" "Uterus" "Vagina" "Whole_blood")

for tissue in ${TISSUES[@]}; do
	mkdir $MAIN_DIR\/$tissue
	n_genes=($(wc -l $EXP_DIR\/$tissue\.rds.erna_exp.norm.sex_age_platform_pcs_peer1_10.txt))
	sbatch \
	--job-name=eRNA\_$tissue \
	--account=g_gamazon_lab \
	--nodes=1 \
	--ntasks=1 \
	--cpus-per-task=2 \
	--mem-per-cpu=300G \
	--time=1-00:00:00 \
	--job-name=$tissue \
	--output $tissue\_erna_predixcan.log  \
	--wrap="
		Rscript $PREDIXCAN_DIR\/predixcan_r.r \
		--model_training \
		--main_dir $MAIN_DIR\/$tissue \
		--plink_file_name $GENO_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.pruned_geno0_maf0.05_hwe1e-6 \
		--expression_file_name $EXP_DIR/$tissue\.rds.erna_exp.norm.sex_age_platform_pcs_peer1_10.txt \
		--n_genes_for_each_subjob $n_genes \
		--annotation_file_name $ANNO_FILE \
		--include_all_gene_types TRUE \
		--parallel
		
		Rscript $PREDIXCAN_DIR\/predixcan_r.r \
         --generate_db_and_cov \
         --main_dir $MAIN_DIR\/$tissue \
         --plink_file_name $GENO_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.pruned_geno0_maf0.05_hwe1e-6 \
         --expression_file_name $EXP_DIR/$tissue\.rds.erna_exp.norm.sex_age_platform_pcs_peer1_10.txt \
         --annotation_file_name $ANNO_FILE \
         --output_file_name $tissue\.erna.predixcan.hg19.exp_cov_norm"
done