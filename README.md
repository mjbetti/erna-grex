# erna-grex [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/mjbetti/erna-grex/blob/master/LICENSE) 
## Introduction  

This repository contains scripts used to run the model training and subsequent analyses described by Betti et al. in "Genetically regulated enhancer RNA expression predicts enhancer-promoter contact frequency and reveals genetic mechanisms at trait-associated loci."

The repository is broken down into the following sub-folders:
* ```grex_model_training```: contains scripts used to train the initial models of genetically regulated enhancer RNA (eRNA) expression (GREx) (**subfigure a**)
* ```contact_frequency_model_training```: contains Jupyter notebooks used to train the deep learning-based models of 3D contact frequency using eRNA and gene GREx in BioVU (**subfigure b**)
* ```trained_contact_models```: contains the optimal deep learning-based models of contact frequency trained in whole blood and cerebellum, respectively, in HDF5 format (**subfigure b**)
* ```schizophrenia```: contains the scripts for running eRNA-based and gene-based TWAS of SCZ and testing the results using MR (**subfigure c and subfigure d**)
* ```uk_biobank```: contains the scripts used to compile the eRNA-based TWAS results across traits in the U.K. Biobank (**subfigure e**)

![alt text](https://github.com/mjbetti/erna-grex/blob/main/Fig1.png?raw=true)

## Citation 
Betti, M.J., Aldrich, M.C., Lin, P., & Gamazon, E.R. (2024). *Genetically regulated enhancer RNA expression predicts enhancer-promoter contact frequency and reveals genetic mechanisms at trait-associated loci*. DOI:  [Preprint](url).

Betti, M.J., Aldrich, M.C., & Gamazon, E.R. (2024). erna-grex datasets (Version 1). Zenodo. DOI: 10.5281/zenodo.11212496.

For questions:  michael.j.betti@vanderbilt.edu, eric.gamazon@vumc.org  

## Data availability
Trained GREx models, as well as all other non-protected data used for model training and other analyses described in the paper can be downloaded from Zenodo at (DOI).

## Dependencies
* For a list of dependencies required to run GREx imputation using the PrediXcan framework, see https://github.com/hakyimlab/MetaXcan.
* For training the deep learning-based contact frequency models and/or running inference using the pre-trained models, the following dependencies were installed:
    * pytorch (1.13.1)
    * torchvision (0.14.1)
    * torchaudio (0.13.1)
    * pytorch-cuda (11.6)
    * cuda-toolkit (11.6)
    * scikit-learn (1.3.1)
    * skorch (0.15.0)
    * pandas (1.5.3)
    * matplotlib (3.7.2)
    * torchmetrics (1.2.0)
    * shap (0.42.1)

Specific versions required may vary, depending on your GPU and which version of CUDA it is using. A .yml file of the full environment used to run these analyses is saved to ```torch-gpu.yml```.

## Usage
### eRNA GREx inference
The following is an example of how one would run an eRNA-based TWAS across 49 tissues:
```
TISSUES=("Adipose_sub" "Adipose_vis" "Adrenal" "Artery_aor" "Artery_cor" "Artery_tib" "Brain_Amy" "Brain_Ant" "Brain_Cau" "Brain_Cerebellum" "Brain_Cer" "Brain_Cortex" "Brain_Fro" "Brain_Hippo" "Brain_Hypo" "Brain_Nuc" "Brain_Put" "Brain_Spi" "Brain_Sub" "Breast" "Cells_ebv" "Cells_fibroblast" "Colon_sigmoid" "Colon_transverse" "Esophagus_gast" "Esophagus_muc" "Esophagus_mus" "Heart_Atr" "Heart_Left" "Kidney" "Liver" "Lung" "Minor_salivary_gland" "Muscle_skeletal" "Nerve" "Ovary" "Pancreas" "Pituitary" "Prostate" "Skin_no_sun" "Skin_sun" "Small_intestine" "Spleen" "Stomach" "Testis" "Thyroid" "Uterus" "Vagina" "Whole_blood")

#Declare the path of the source directory containing the MetaXcan scripts, as well as the directories of the inputs, outputs, and GTEx models
SCRIPT_DIR=/data/aldrich_lab/bin/MetaXcan_scripts/MetaXcan/software
MAIN_DIR=/data/g_gamazon_lab/bettimj/eRNA_predixcan_models/attempt4/training/all_models
GWAS_DIR=/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/for_metaxcan_in
OUT_DIR=/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/schizophrenia/twas_all

#Run MetaXcan using the prepared GTEx v8 models and AALC GWAS results
for tissue in ${TISSUES[@]}; do
	python $SCRIPT_DIR\/SPrediXcan.py \
	--model_db_path $MAIN_DIR\/$tissue\.erna.predixcan.hg19.exp_cov_norm.db \
	--covariance $MAIN_DIR\/$tissue\.erna.predixcan.hg19.exp_cov_norm.cov \
	--gwas_folder $GWAS_DIR \
	--snp_column SNP \
	--effect_allele_column A2 \
	--non_effect_allele_column A1 \
	--pvalue_column PVAL \
	--beta_column BETA \
	--output_file $OUT_DIR\/schizophrenia2022_eRNA_imputations.$tissue\.csv \
	--keep_non_rsid
done
```

### Enhancer-gene contact frequency prediction
The following is an example of how to predict pairwise contact frequencies using input eRNA and gene GREx values. If using one of the pre-trained models, we recommend using the model trained on cerebellum, as it shows robust cross-tissue portability compared with the model trained on whole blood.
```
import torch
import shap

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

class Net(nn.Module):
    def __init__(self, input_size=2, hidden_size=90, output_size=1, num_hidden_layers=2, weight_init_hidden=init.xavier_normal_, weight_init_out=init.xavier_uniform_):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.num_hidden_layers = num_hidden_layers - 1
        self.hidden_layers = nn.ModuleList([nn.Linear(hidden_size, hidden_size) for i in range(num_hidden_layers - 1)])
        self.fc2 = nn.Linear(hidden_size, output_size)
        #self.activation = activation()
        # manually init weights
        weight_init_hidden(self.fc1.weight)
        for i, hidden_layer in enumerate(self.hidden_layers):
            weight_init_hidden(hidden_layer.weight)
        weight_init_out(self.fc2.weight)

    def forward(self, x):
        x = torch.nn.functional.softsign(self.fc1(x))
        for i in range(self.num_hidden_layers):
            x = torch.nn.functional.softsign(self.hidden_layers[i](x))
        x = self.fc2(x)
        return x
    
model = Net().to(device)
model.load_state_dict(torch.load('cerebellum_best_contact_regression_model.h5'))

y_pred = best(test_set_final.to("cuda:0"))
```
