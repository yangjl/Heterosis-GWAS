# Jinliang Yang
# Purpose: control GenSel running
# date: June/20/2014

source("~/Documents/Heterosis_GWAS/HGWAS_proj/lib/sh4GenSel.R")
### real runs:
sh4GenSel(pwd="~/Documents/Heterosis_GWAS/Method/GenSel", 
          sh="KC_run41000.sh",
          geno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/geno_kc_1.3M_GenSel.newbin", 
          pheno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/kc_v2_GenSel_fullset.txt",
          chainLength=41000, burnin=1000, varGenotypic=2795, varResidual=4008)

sh4GenSel(pwd="~/Documents/Heterosis_GWAS/Method/GenSel", 
          sh="CL_run41000.sh",
          geno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/geno_cl_1.3M_GenSel.newbin", 
          pheno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/cl_v2_GenSel_fullset.txt",
          chainLength=41000, burnin=1000, varGenotypic=190, varResidual=143)


sh4GenSel(pwd="~/Documents/Heterosis_GWAS/Method/GenSel", 
          sh="CD_run41000.sh",
          geno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/geno_cd_1.3M_GenSel.newbin", 
          pheno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/cd_v2_GenSel_fullset.txt",
          chainLength=41000, burnin=1000, varGenotypic=5, varResidual=2)

sh4GenSel(pwd="~/Documents/Heterosis_GWAS/Method/GenSel", 
          sh="CW_run41000.sh",
          geno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/geno_cw_1.3M_GenSel.newbin", 
          pheno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/cw_v2_GenSel_fullset.txt",
          chainLength=41000, burnin=1000, varGenotypic=21, varResidual=11)


sh4GenSel(pwd="~/Documents/Heterosis_GWAS/Method/GenSel", 
          sh="AKW_run41000.sh",
          geno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/geno_akw_1.3M_GenSel.newbin", 
          pheno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/akw_v2_GenSel_fullset.txt",
          chainLength=41000, burnin=1000, varGenotypic=0.0004, varResidual=0.004)

sh4GenSel(pwd="~/Documents/Heterosis_GWAS/Method/GenSel", 
          sh="TKW_run41000.sh",
          geno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/geno_tkw_1.3M_GenSel.newbin", 
          pheno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/tkw_v2_GenSel_fullset.txt",
          chainLength=41000, burnin=1000, varGenotypic=294, varResidual=282)

sh4GenSel(pwd="~/Documents/Heterosis_GWAS/Method/GenSel", 
          sh="KRN_run41000.sh",
          geno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/geno_cd_1.3M_GenSel.newbin", 
          pheno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/krn_v2_GenSel_fullset.txt",
          chainLength=41000, burnin=1000, varGenotypic=2.1, varResidual=1)

####stop here!!!
####DONE