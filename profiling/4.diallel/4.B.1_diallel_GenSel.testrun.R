# Jinliang Yang
# Purpose: control GenSel running
# date: June/20/2014

source("~/Documents/Heterosis_GWAS/HGWAS_proj/lib/sh4GenSel.R")
### real runs:
pwd = "/home/NSF-SAM-GS/HGWAS/SNPdiallel/dtestrun/"
traits <- c("krn", "cd", "cl", "cw", "akw", "tkw", "kc")
genovar <- c(2.1, 5, 190, 21, 0.0004, 294, 2795)
residual <- c(1, 2, 143, 11, 0.004, 282, 4008)

################ testrun for add genotype
for(i in 1:7){
  mysh <- paste(traits[i], "_add_run1000.sh", sep="")
  mygeno <- paste("/home/NSF-SAM-GS/HGWAS/SNPdiallel/snp", traits[i], "_1m_add_gensel", sep="")
  mypheno <- paste("/mnt/02/yangjl/Documents/Heterosis_GWAS/HGWAS_proj/data/", traits[i],
                   "_v2_GenSel_diallel.txt", sep="")
  sh4GenSel(pwd=pwd, pai=0.999, chainLength=1000, burnin=100, 
            varGenotypic=genovar[i], varResidual= residual[i],
            sh=mysh, geno=mygeno, pheno=mypheno) 
}

################ testrun for dom genotype
for(i in 1:7){
  mysh <- paste(traits[i], "_dom_run1000.sh", sep="")
  mygeno <- paste("/home/NSF-SAM-GS/HGWAS/SNPdiallel/snp", traits[i], "_1m_dom_gensel", sep="")
  mypheno <- paste("/mnt/02/yangjl/Documents/Heterosis_GWAS/HGWAS_proj/data/", traits[i],
                   "_v2_GenSel_diallel.txt", sep="")
  sh4GenSel(pwd=pwd, pai=0.9995, chainLength=1000, burnin=100, 
            varGenotypic=genovar[i], varResidual= residual[i],
            sh=mysh, geno=mygeno, pheno=mypheno) 
}

