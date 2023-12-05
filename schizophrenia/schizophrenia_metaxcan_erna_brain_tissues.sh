TISSUES=("Brain_Amy" "Brain_Ant" "Brain_Cau" "Brain_Cerebellum" "Brain_Cer" "Brain_Cortex" "Brain_Fro" "Brain_Hippo" "Brain_Hypo" "Brain_Nuc" "Brain_Put" "Brain_Spi" "Brain_Sub")

#Declare the path of the source directory containing the MetaXcan scripts, as well as the directories of the inputs, outputs, and GTEx models
SCRIPT_DIR=/data/aldrich_lab/bin/MetaXcan_scripts/MetaXcan/software
MAIN_DIR=/home/bettimj/gamazon_rotation/eRNA_predixcan_models
GWAS_DIR=/home/bettimj/gamazon_rotation/eRNA_predixcan_models/gwas_sum_stats/schizophrenia/liftover_hg38/for_metaxcan_in
OUT_DIR=/home/bettimj/gamazon_rotation/eRNA_predixcan_models/gwas_sum_stats/schizophrenia/imputations

#Run MetaXcan using the prepared GTEx v8 models and AALC GWAS results
for tissue in ${TISSUES[@]}; do
	python $SCRIPT_DIR\/SPrediXcan.py \
	--model_db_path $MAIN_DIR\/$tissue\/output/$tissue\.eRNA.db \
	--covariance $MAIN_DIR\/$tissue\/output/$tissue\.eRNA.cov \
	--gwas_folder $GWAS_DIR \
	--snp_column varID \
	--effect_allele_column A1 \
	--non_effect_allele_column A2 \
	--pvalue_column P \
	--or_column OR \
	--output_file $OUT_DIR\/schizophrenia_eRNA_imputations.$tissue\.csv
done
