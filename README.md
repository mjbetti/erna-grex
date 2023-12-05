# erna-grex [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/mjbetti/erna-grex/blob/master/LICENSE) 
## Introduction  

This repository contains scripts used to run the model training and subsequent analyses described by Betti et al. in "Genetically regulated enhancer RNA expression is predictive of enhancer-gene contact frequency and reveals genetic mechanisms at trait-associated loci."

The repository is broken down into the following sub-folders:
* ```grex_model_training```: contains scripts used to train the initial models of genetically regulated enhancer RNA (eRNA) expression (GREx)
* ```contact_frequency_model_training```: contains Jupyter notebooks used to train the deep learning-based models of 3D contact frequency using eRNA and gene GREx in BioVU
* ```trained_contact_models```: contains the optimal deep learning-based models of contact frequency trained in whole blood and cerebellum, respectively, in HDF5 format
* ```schizophrenia```: contains the scripts for running eRNA-based and gene-based TWAS of SCZ and testing the results using MR
* ```uk_biobank```: contains the scripts used to compile the eRNA-based TWAS results across traits in the U.K. Biobank

![alt text](https://github.com/mjbetti/erna-grex/blob/main/FigAbs.png?raw=true)

## Citation 
Betti, M.J., Aldrich, M.C., Lin, P., & Gamazon, E.R. (2023). *Genetically regulated enhancer RNA expression is predictive of enhancer-gene contact frequency and reveals genetic mechanisms at trait-associated loci*. doi:  [Preprint](url).

Betti, M.J., Aldrich, M.C., & Gamazon, E.R. (2023). erna-grex datasets (Version 1). Zenodo. DOI

For questions:  michael.j.betti@vanderbilt.edu, eric.gamazon@vumc.org  

## Data availability
