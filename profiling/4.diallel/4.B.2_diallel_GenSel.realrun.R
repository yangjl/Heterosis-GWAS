# Jinliang Yang
# Purpose: control GenSel running
# date: June/20/2014

source("~/Documents/Heterosis_GWAS/HGWAS_proj/lib/sh4GenSel.R")
### real runs:
traits <- c("krn", "cd", "cl", "cw", "akw", "tkw", "kc")
pwd = "/home/NSF-SAM-GS/HGWAS/SNPdiallel/drealrun/"

################ testrun for dom genotype
genovar <- c(1.9, 7.6, 373, 91, 0.0017, 2515, 15485)
residual <- c(0.3, 0.7, 105, 7, 0.0007, 178, 2258)

for(i in 1:7){
  mysh <- paste(traits[i], "_dom_run41000.sh", sep="")
  mygeno <- paste("/home/NSF-SAM-GS/HGWAS/SNPdiallel/snp", traits[i], 
                  "_1m_dom_gensel.newbin", sep="")
  mypheno <- paste("/mnt/02/yangjl/Documents/Heterosis_GWAS/HGWAS_proj/data/", traits[i],
                   "_v2_GenSel_diallel.txt", sep="")
  sh4GenSel(pwd=pwd, pai=0.9995, chainLength=41000, burnin=1000, 
            varGenotypic=genovar[i], varResidual= residual[i],
            sh=mysh, geno=mygeno, pheno=mypheno) 
}


################ testrun for add genotype
traits <- c("krn", "cd", "cl", "cw", "akw", "tkw", "kc")
pwd = "/home/NSF-SAM-GS/HGWAS/SNPdiallel/drealrun/"
genovar <- c(1.7, 6.5, 198, 70, 0.0012, 1474, 9673)
residual <- c(0.5, 1.6, 280, 26, 0.0012, 1182, 7901)

for(i in 1:7){
  mysh <- paste(traits[i], "_add_run41000.sh", sep="")
  mygeno <- paste("/home/NSF-SAM-GS/HGWAS/SNPdiallel/snp", traits[i], 
                  "_1m_add_gensel.newbin", sep="")
  mypheno <- paste("/mnt/02/yangjl/Documents/Heterosis_GWAS/HGWAS_proj/data/", traits[i],
                   "_v2_GenSel_diallel.txt", sep="")
  sh4GenSel(pwd=pwd, pai=0.999, chainLength=41000, burnin=1000, 
            varGenotypic=genovar[i], varResidual= residual[i],
            sh=mysh, geno=mygeno, pheno=mypheno) 
}


