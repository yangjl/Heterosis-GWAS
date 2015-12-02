# Jinliang Yang
# Purpose: Prediction of SAM large and SAM small
# date: July.14.2014
# location: server.9



setupPredict <- function(trait=traits[1], pwd=pwd, mysh="prediction_s2.sh"){
  for(j in 1:1){ ### change j to add more randomization 
    
    #### add 
    for(k in 1:5){
      myinp <- paste(trait, "_", j, "_", k, "_pred_add.inp", sep="")
      mypheno1 <- paste("/home/NSF-SAM-GS/HGWAS/SNPdiallel/dcv2/pheno/", trait,
                        "_", j, "_val_s", k, ".txt", sep="")
      mygeno1 <- paste("/home/NSF-SAM-GS/HGWAS/SNPdiallel/snp", trait, 
                       "_1m_add_gensel.newbin", sep="")
      #map="/home/NSF-SAM-GS/snpmatrix/SAM_SNPs_chrall_samid.map"
      mymrkRes <- paste("/home/NSF-SAM-GS/HGWAS/SNPdiallel/dcv2/", trait,
                        "/", trait, "_", j, "_", k, "_train_add.mrkRes1", sep="")
      
      sh4Predict_v2(pwd=pwd, sh=mysh, inp=myinp,
                 geno=mygeno1, pheno=mypheno1, mrkRes=mymrkRes, shappend=TRUE) 
      
    }
    
    ####### dom
    for(k in 1:5){
      myinp2 <- paste(trait, "_", j, "_", k, "_pred_dom.inp", sep="")
      mypheno2 <- paste("/home/NSF-SAM-GS/HGWAS/SNPdiallel/dcv2/pheno/", trait,
                        "_", j, "_val_s", k, ".txt", sep="")
      mygeno2 <- paste("/home/NSF-SAM-GS/HGWAS/SNPdiallel/snp", trait, 
                       "_1m_dom_gensel.newbin", sep="")
      #map="/home/NSF-SAM-GS/snpmatrix/SAM_SNPs_chrall_samid.map"
      mymrkRes2 <- paste("/home/NSF-SAM-GS/HGWAS/SNPdiallel/dcv2/", trait,
                        "/", trait, "_", j, "_", k, "_train_dom.mrkRes1", sep="")
      sh4Predict_v2(pwd=pwd, sh=mysh, inp=myinp2,
                 geno=mygeno2, pheno=mypheno2, mrkRes=mymrkRes2, shappend=TRUE) 
    }
  }
}
################
source("~/Documents/Heterosis_GWAS/HGWAS_proj/lib/sh4Predict_v2.R")

pwd = "/home/NSF-SAM-GS/HGWAS/SNPdiallel/predict/s2_predict"

### real runs:
traits <- c("krn", "cd", "cl", "cw", "akw", "tkw", "kc")
setupPredict(trait=traits[1], pwd=pwd, mysh="prediction_s2.sh")
setupPredict(trait=traits[2], pwd=pwd, mysh="prediction_s2.sh")
setupPredict(trait=traits[3], pwd=pwd, mysh="prediction_s2.sh")
setupPredict(trait=traits[4], pwd=pwd, mysh="prediction_s2.sh")
setupPredict(trait=traits[5], pwd=pwd, mysh="prediction_s2.sh")
setupPredict(trait=traits[6], pwd=pwd, mysh="prediction_s2.sh")
setupPredict(trait=traits[7], pwd=pwd, mysh="prediction_s2.sh")
