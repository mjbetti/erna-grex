#Declare the paths of the VCF and pop files, as well as the directory containing the JTI SQL databases
VCF=/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training/afr_cov/cov_matrices/recode_oddball/SCCS_gtexIDs_v3.withchr.hg19.vcf
POP=/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training/afr_cov/recoded_bim/samples.txt
DB_DIR=/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training/all_models

TISSUES=("Adipose_sub" "Adipose_vis" "Adrenal" "Artery_aor" "Artery_cor" "Artery_tib" "Brain_Amy" "Brain_Ant" "Brain_Cau" "Brain_Cerebellum" "Brain_Cer" "Brain_Cortex" "Brain_Fro" "Brain_Hippo" "Brain_Hypo" "Brain_Nuc" "Brain_Put" "Brain_Spi" "Brain_Sub" "Breast" "Cells_ebv" "Cells_fibroblast" "Colon_sigmoid" "Colon_transverse" "Esophagus_gast" "Esophagus_muc" "Esophagus_mus" "Heart_Atr" "Heart_Left" "Kidney" "Liver" "Lung" "Minor_salivary_gland" "Muscle_skeletal" "Nerve" "Ovary" "Pancreas" "Pituitary" "Prostate" "Skin_no_sun" "Skin_sun" "Small_intestine" "Spleen" "Stomach" "Testis" "Thyroid" "Uterus" "Vagina" "Whole_blood")

#Loop through each tissue, submitting a SLURM job for each

for TISSUE in "${TISSUES[@]}"; do
	sbatch \
	--nodes=1 \
	--ntasks=1 \
	--cpus-per-task=1 \
	--mem-per-cpu=64G \
	--time=5-00:00:00 \
	--job-name=$TISSUE\_mxcovar \
	--account aldrich_lab \
	--wrap="python /home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training/afr_cov/cov_matrices/mxcovar_non1KG_erna.py \
	--vcf /home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training/afr_cov/cov_matrices/recode_oddball/SCCS_gtexIDs_v3.withchr.hg19.vcf \
	--pop /home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training/afr_cov/recoded_bim/samples.txt \
	$DB_DIR\/$TISSUE\.erna.predixcan.hg19.exp_cov_norm.db"
done