#Declare the paths of the liftOver directory, as well as the working directory
LIFTOVER_BIN=/home/bettimj/liftOver
CHAIN=/home/bettimj/reference_genomes/hg19ToHg38.over.chain.gz
WORKING_DIR=/home/bettimj/gamazon_rotation/eRNA_predixcan_models/gwas_sum_stats/schizophrenia/liftover_hg38

#Lift over the bed coordinates from hg19 to hg38
$LIFTOVER_BIN \
$WORKING_DIR\/schizophrenia_daner_PGC_SCZ52_0513a.hg19.bed \
$CHAIN \
$WORKING_DIR\/schizophrenia_daner_PGC_SCZ52_0513a.hg38.bed \
$WORKING_DIR\/schizophrenia_daner_PGC_SCZ52_0513a.hg19.unlifted.bed