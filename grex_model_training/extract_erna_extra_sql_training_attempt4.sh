#Declare the path of the root directory, as well as the desired output directory
ROOT_DIR=/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt4/training
OUT_DIR=$ROOT_DIR\/tissue_extra_r2_lists

#Create an array containing the names of all tissue types for which eRNA models were trained
TISSUES=("Adipose_sub" "Adipose_vis" "Adrenal" "Artery_aor" "Artery_cor" "Artery_tib" "Brain_Amy" "Brain_Ant" "Brain_Cau" "Brain_Cerebellum" "Brain_Cer" "Brain_Cortex" "Brain_Fro" "Brain_Hippo" "Brain_Hypo" "Brain_Nuc" "Brain_Put" "Brain_Spi" "Brain_Sub" "Cells_ebv" "Cells_fibroblast" "Cells_leukemia" "Colon_sigmoid" "Colon_transverse" "Esophagus_gast" "Esophagus_muc" "Esophagus_mus" "Heart_Atr" "Heart_Left" "Kidney" "Liver" "Lung" "Minor_salivary_gland" "Muscle_skeletal" "Nerve" "Ovary" "Pancreas" "Pituitary" "Prostate" "Skin_no_sun" "Skin_sun" "Small_intestine" "Spleen" "Stomach" "Testis" "Thyroid" "Uterus" "Vagina" "Whole_blood")

# for TISSUE in ${TISSUES[@]}; do
# 	sqlite3 $ROOT_DIR\/$TISSUE\/output/$TISSUE\.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
# 	.output $tissue.erna.extra.predixcan.hg19.txt
# 	select * from extra;
# 	END_SQL
# done

#Adipose_sub
sqlite3 $ROOT_DIR\/Adipose_sub\/output/Adipose_sub.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Adipose_sub.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Adipose_vis
sqlite3 $ROOT_DIR\/Adipose_vis\/output/Adipose_vis.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Adipose_vis.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Adrenal
sqlite3 $ROOT_DIR\/Adrenal\/output/Adrenal.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Adrenal.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Artery_aor
sqlite3 $ROOT_DIR\/Artery_aor\/output/Artery_aor.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Artery_aor.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Artery_cor
sqlite3 $ROOT_DIR\/Artery_cor\/output/Artery_cor.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Artery_cor.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Artery_tib
sqlite3 $ROOT_DIR\/Artery_tib\/output/Artery_tib.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Artery_tib.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Bladder
sqlite3 $ROOT_DIR\/Bladder\/output/Bladder.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Bladder.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Brain_Amy
sqlite3 $ROOT_DIR\/Brain_Amy\/output/Brain_Amy.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Brain_Amy.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Brain_Ant
sqlite3 $ROOT_DIR\/Brain_Ant\/output/Brain_Ant.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Brain_Ant.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Brain_Cau
sqlite3 $ROOT_DIR\/Brain_Cau\/output/Brain_Cau.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Brain_Cau.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Brain_Cer
sqlite3 $ROOT_DIR\/Brain_Cer\/output/Brain_Cer.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Brain_Cer.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Brain_Cerebellum
sqlite3 $ROOT_DIR\/Brain_Cerebellum\/output/Brain_Cerebellum.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Brain_Cerebellum.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Brain_Cortex
sqlite3 $ROOT_DIR\/Brain_Cortex\/output/Brain_Cortex.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Brain_Cortex.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Brain_Fro
sqlite3 $ROOT_DIR\/Brain_Fro\/output/Brain_Fro.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Brain_Fro.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Brain_Hippo
sqlite3 $ROOT_DIR\/Brain_Hippo\/output/Brain_Hippo.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Brain_Hippo.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Brain_Hypo
sqlite3 $ROOT_DIR\/Brain_Hypo\/output/Brain_Hypo.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Brain_Hypo.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Brain_Nuc
sqlite3 $ROOT_DIR\/Brain_Nuc\/output/Brain_Nuc.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Brain_Nuc.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Brain_Put
sqlite3 $ROOT_DIR\/Brain_Put\/output/Brain_Put.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Brain_Put.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Brain_Spi
sqlite3 $ROOT_DIR\/Brain_Spi\/output/Brain_Spi.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Brain_Spi.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Brain_Sub
sqlite3 $ROOT_DIR\/Brain_Sub\/output/Brain_Sub.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Brain_Sub.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Breast
sqlite3 $ROOT_DIR\/Breast\/output/Breast.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Breast.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Cells_ebv
sqlite3 $ROOT_DIR\/Cells_ebv\/output/Cells_ebv.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Cells_ebv.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Cells_fibroblast
sqlite3 $ROOT_DIR\/Cells_fibroblast\/output/Cells_fibroblast.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Cells_fibroblast.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Cells_leukemia
# sqlite3 $ROOT_DIR\/Cells_leukemia\/output/Cells_leukemia.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
# .output Cells_leukemia.erna.extra.predixcan.hg19.txt
# select * from extra;
# END_SQL

#Cervix_ecto
sqlite3 $ROOT_DIR\/Cervix_ecto\/output/Cervix_ecto.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Cervix_ecto.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Cervix_endo
sqlite3 $ROOT_DIR\/Cervix_endo\/output/Cervix_endo.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Cervix_endo.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Colon_sigmoid
sqlite3 $ROOT_DIR\/Colon_sigmoid\/output/Colon_sigmoid.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Colon_sigmoid.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Colon_transverse
sqlite3 $ROOT_DIR\/Colon_transverse\/output/Colon_transverse.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Colon_transverse.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Esophagus_gast
sqlite3 $ROOT_DIR\/Esophagus_gast\/output/Esophagus_gast.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Esophagus_gast.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Esophagus_muc
sqlite3 $ROOT_DIR\/Esophagus_muc\/output/Esophagus_muc.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Esophagus_muc.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Esophagus_mus
sqlite3 $ROOT_DIR\/Esophagus_mus\/output/Esophagus_mus.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Esophagus_mus.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Fallopian_tube
sqlite3 $ROOT_DIR\/Fallopian_tube\/output/Fallopian_tube.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Fallopian_tube.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Heart_Atr
sqlite3 $ROOT_DIR\/Heart_Atr\/output/Heart_Atr.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Heart_Atr.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Heart_Left
sqlite3 $ROOT_DIR\/Heart_Left\/output/Heart_Left.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Heart_Left.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Kidney
sqlite3 $ROOT_DIR\/Kidney\/output/Kidney.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Kidney.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Liver
sqlite3 $ROOT_DIR\/Liver\/output/Liver.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Liver.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Lung
sqlite3 $ROOT_DIR\/Lung\/output/Lung.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Lung.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Minor_salivary_gland
sqlite3 $ROOT_DIR\/Minor_salivary_gland\/output/Minor_salivary_gland.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Minor_salivary_gland.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Muscle_skeletal
sqlite3 $ROOT_DIR\/Muscle_skeletal\/output/Muscle_skeletal.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Muscle_skeletal.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Nerve
sqlite3 $ROOT_DIR\/Nerve\/output/Nerve.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Nerve.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Ovary
sqlite3 $ROOT_DIR\/Ovary\/output/Ovary.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Ovary.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Pancreas
sqlite3 $ROOT_DIR\/Pancreas\/output/Pancreas.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Pancreas.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Pituitary
sqlite3 $ROOT_DIR\/Pituitary\/output/Pituitary.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Pituitary.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Prostate
sqlite3 $ROOT_DIR\/Prostate\/output/Prostate.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Prostate.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Skin_no_sun
sqlite3 $ROOT_DIR\/Skin_no_sun\/output/Skin_no_sun.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Skin_no_sun.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Skin_sun
sqlite3 $ROOT_DIR\/Skin_sun\/output/Skin_sun.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Skin_sun.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Small_intestine
sqlite3 $ROOT_DIR\/Small_intestine\/output/Small_intestine.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Small_intestine.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Spleen
sqlite3 $ROOT_DIR\/Spleen\/output/Spleen.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Spleen.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Stomach
sqlite3 $ROOT_DIR\/Stomach\/output/Stomach.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Stomach.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Testis
sqlite3 $ROOT_DIR\/Testis\/output/Testis.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Testis.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Thyroid
sqlite3 $ROOT_DIR\/Thyroid\/output/Thyroid.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Thyroid.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Uterus
sqlite3 $ROOT_DIR\/Uterus\/output/Uterus.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Uterus.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Vagina
sqlite3 $ROOT_DIR\/Vagina\/output/Vagina.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Vagina.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL

#Whole_blood
sqlite3 $ROOT_DIR\/Whole_blood\/output/Whole_blood.erna.predixcan.hg19.exp_cov_norm.db << 'END_SQL'
.output Whole_blood.erna.extra.predixcan.hg19.txt
select * from extra;
END_SQL