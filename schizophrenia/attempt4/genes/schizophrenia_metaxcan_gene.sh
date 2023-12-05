TISSUES=("Adipose_Subcutaneous" "Adipose_Visceral_Omentum" "Adrenal_Gland" "Artery_Aorta" "Artery_Coronary" "Artery_Tibial" "Brain_Amygdala" "Brain_Anterior_cingulate_cortex_BA24" "Brain_Caudate_basal_ganglia" "Brain_Cerebellar_Hemisphere" "Brain_Cerebellum" "Brain_Cortex" "Brain_Frontal_Cortex_BA9" "Brain_Hippocampus" "Brain_Hypothalamus" "Brain_Nucleus_accumbens_basal_ganglia" "Brain_Putamen_basal_ganglia" "Brain_Spinal_cord_cervical_c-1" "Brain_Substantia_nigra" "Breast_Mammary_Tissue" "Cells_Cultured_fibroblasts" "Cells_EBV-transformed_lymphocytes" "Colon_Sigmoid" "Colon_Transverse" "Esophagus_Gastroesophageal_Junction" "Esophagus_Mucosa" "Esophagus_Muscularis" "Heart_Atrial_Appendage" "Heart_Left_Ventricle" "Kidney_Cortex" "Liver" "Lung" "Minor_Salivary_Gland" "Muscle_Skeletal" "Nerve_Tibial" "Ovary" "Pancreas" "Pituitary" "Prostate" "Skin_Not_Sun_Exposed_Suprapubic" "Skin_Sun_Exposed_Lower_leg" "Small_Intestine_Terminal_Ileum" "Spleen" "Stomach" "Testis" "Thyroid" "Uterus" "Vagina" "Whole_Blood"
)
#Declare the path of the source directory containing the MetaXcan scripts, as well as the directories of the inputs, outputs, and GTEx models
SCRIPT_DIR=/data/aldrich_lab/bin/MetaXcan_scripts/MetaXcan/software
MAIN_DIR=/home/bettimj/gamazon_rotation/gtex8_models/predixcan
GWAS_DIR=/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/for_metaxcan_in
OUT_DIR=/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all_genes

#Run MetaXcan using the prepared models and GWAS results
for tissue in ${TISSUES[@]}; do
python $SCRIPT_DIR\/SPrediXcan.py \
--model_db_path $MAIN_DIR\/PrediXcan_$tissue\.db \
--covariance $MAIN_DIR\/PrediXcan_$tissue\.txt.gz \
--gwas_folder $GWAS_DIR \
--snp_column ID \
--effect_allele_column A2 \
--non_effect_allele_column A1 \
--pvalue_column PVAL \
--beta_column BETA \
--output_file $OUT_DIR\/schizophrenia2022_gene_imputations.$tissue\.csv \
--keep_non_rsid
done