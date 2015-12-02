# Jinliang Yang
# Purpose: cross-validation
# date: 9/13/2014

setupGenSel_s1 <- function(trait=traits[1], i=1){
  for(j in 1:1){ ### change j to add more randomization 
        
    mysh <- paste(trait, j, "run41000_training_s1.sh", sep="_")
    
    #### add 
    for(k in 1:5){
      mypheno1 <- paste("/home/NSF-SAM-GS/HGWAS/SNPdiallel/dcv1/pheno/", trait,
                        "_", j, "_train", k, ".txt", sep="")
      mygeno1 <- paste("/home/NSF-SAM-GS/HGWAS/SNPdiallel/snp", trait, 
                       "_1m_dom_gensel.newbin", sep="")
      inpfile <- paste(trait, j, k, "train_dom.inp", sep="_")
      sh4GenSel(pwd=pwd, inp=inpfile,
                pai=0.9995, chainLength=11000, burnin=1000, 
                varGenotypic=genovar[i], varResidual= residual[i],
                sh=mysh, geno=mygeno1, pheno=mypheno1, shappend=TRUE) 
    }
    
    ####### dom
    for(k in 1:5){
      mypheno2 <- paste("/home/NSF-SAM-GS/HGWAS/SNPdiallel/dcv1/pheno/", trait,
                        "_", j, "_train", k, ".txt", sep="")
      mygeno2 <- paste("/home/NSF-SAM-GS/HGWAS/SNPdiallel/snp", trait, 
                       "_1m_add_gensel.newbin", sep="")
      
      inpfile <- paste(trait, j, k, "train_add.inp", sep="_")
      sh4GenSel(pwd=pwd, inp=inpfile,
                pai=0.999, chainLength=11000, burnin=1000, 
                varGenotypic=genovar[i], varResidual= residual[i],
                sh=mysh, geno=mygeno2, pheno=mypheno2, shappend=TRUE) 
    }
  }
}
################
source("~/Documents/Heterosis_GWAS/HGWAS_proj/lib/sh4GenSelv2.R")
### real runs:
traits <- c("krn", "cd", "cl", "cw", "akw", "tkw", "kc")


################ testrun for dom genotype
genovar <- c(1.9, 7.6, 373, 91, 0.0017, 2515, 15485)
residual <- c(0.3, 0.7, 105, 7, 0.0007, 178, 2258)

pwd = "/home/NSF-SAM-GS/HGWAS/SNPdiallel/dcv1/krn"
setupGenSel_s1(trait=traits[1], i=1)
pwd = "/home/NSF-SAM-GS/HGWAS/SNPdiallel/dcv1/cd"
setupGenSel_s1(trait=traits[2], i=2)
pwd = "/home/NSF-SAM-GS/HGWAS/SNPdiallel/dcv1/cl"
setupGenSel_s1(trait=traits[3], i=3)
pwd = "/home/NSF-SAM-GS/HGWAS/SNPdiallel/dcv1/cw"
setupGenSel_s1(trait=traits[4], i=4)
pwd = "/home/NSF-SAM-GS/HGWAS/SNPdiallel/dcv1/akw"
setupGenSel_s1(trait=traits[5], i=5)
pwd = "/home/NSF-SAM-GS/HGWAS/SNPdiallel/dcv1/tkw"
setupGenSel_s1(trait=traits[6], i=6)
pwd = "/home/NSF-SAM-GS/HGWAS/SNPdiallel/dcv1/kc"
setupGenSel_s1(trait=traits[7], i=7)


