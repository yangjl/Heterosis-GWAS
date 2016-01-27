# Jinliang Yang
# Purpose: control GenSel running
# date: June/20/2014

source("~/Documents/Heterosis_GWAS/HGWAS_proj/lib/sh4GenSel.R")
### test runs:
sh4GenSel(pwd="~/Documents/Heterosis_GWAS/Method/GenSel/testrun", 
          sh="KC_test.sh",
          geno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/geno_kc_1.3M_GenSel.newbin", 
          pheno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/kc_v2_GenSel_fullset.txt",
          chainLength=1000, burnin=100, varGenotypic=1, varResidual=2)

sh4GenSel(pwd="~/Documents/Heterosis_GWAS/Method/GenSel/testrun", 
          sh="CL_test.sh",
          geno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/geno_cl_1.3M_GenSel.newbin", 
          pheno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/cl_v2_GenSel_fullset.txt",
          chainLength=1000, burnin=100, varGenotypic=1, varResidual=2)

sh4GenSel(pwd="~/Documents/Heterosis_GWAS/Method/GenSel/testrun", 
          sh="CD_test.sh",
          geno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/geno_cd_1.3M_GenSel.newbin", 
          pheno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/cd_v2_GenSel_fullset.txt",
          chainLength=1000, burnin=100, varGenotypic=1, varResidual=2)

sh4GenSel(pwd="~/Documents/Heterosis_GWAS/Method/GenSel/testrun", 
          sh="CW_test.sh",
          geno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/geno_cw_1.3M_GenSel.newbin", 
          pheno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/cw_v2_GenSel_fullset.txt",
          chainLength=1000, burnin=100, varGenotypic=1, varResidual=2)

sh4GenSel(pwd="~/Documents/Heterosis_GWAS/Method/GenSel/testrun", 
          sh="AKW_test.sh",
          geno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/geno_akw_1.3M_GenSel.newbin", 
          pheno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/akw_v2_GenSel_fullset.txt",
          chainLength=1000, burnin=100, varGenotypic=1, varResidual=2)

sh4GenSel(pwd="~/Documents/Heterosis_GWAS/Method/GenSel/testrun", 
          sh="TKW_test.sh",
          geno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/geno_tkw_1.3M_GenSel.newbin", 
          pheno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/tkw_v2_GenSel_fullset.txt",
          chainLength=1000, burnin=100, varGenotypic=1, varResidual=2)

sh4GenSel(pwd="~/Documents/Heterosis_GWAS/Method/GenSel/testrun", 
          sh="KRN_test.sh",
          geno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/geno_cd_1.3M_GenSel.newbin", 
          pheno="/mnt/02/yangjl/Documents/Heterosis_GWAS/Method/Stepwise/krn_v2_GenSel_fullset.txt",
          chainLength=1000, burnin=100, varGenotypic=1, varResidual=2)
