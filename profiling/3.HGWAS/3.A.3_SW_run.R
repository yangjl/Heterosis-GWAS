# Jinliang Yang
# Purpose: control GenSel running
# date: June/20/2014

source("~/Documents/Heterosis_GWAS/HGWAS_proj/lib/sh4Stepwise.R")
### runs:


sh4Stepwise(pwd="~/Documents/Heterosis_GWAS/Method/Stepwise/", 
            sh="KRN_sw_v3.sh", email=TRUE,
            geno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/geno_cd_1.3M_GenSel.newbin", 
            pheno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/krn_v2_GenSel_fullset.txt"
            )

####stop here!!!
####DONE