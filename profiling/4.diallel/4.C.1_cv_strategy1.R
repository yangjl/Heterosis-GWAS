# Jinliang Yang
# 5 fold traits split
# date: June/20/2014

### trait split
traits <- c("krn", "cd", "cl", "cw", "akw", "tkw", "kc")
phenopwd <- "/home/NSF-SAM-GS/HGWAS/SNPdiallel/dcv1/pheno"

splitPheno(traits, splittimes=20, outputpwd=phenopwd)


############################  
splitPheno <- function(traits, splittimes=20, outputpwd="~"){
  
  mypheno <- paste("/mnt/02/yangjl/Documents/Heterosis_GWAS/HGWAS_proj/data/", "krn",
                   "_v2_GenSel_diallel.txt", sep="")
  pheno <- read.table(mypheno, header=TRUE)
  
  
  setwd(outputpwd)
  for(i in 1:splittimes){
    idx <- sample(1:nrow(pheno), nrow(pheno))
    myf <- round(nrow(pheno)/5, 0)
    
    for(j in 1:length(traits)){
      mypheno <- paste("/mnt/02/yangjl/Documents/Heterosis_GWAS/HGWAS_proj/data/", traits[j],
                       "_v2_GenSel_diallel.txt", sep="")
      pheno <- read.table(mypheno, header=TRUE)
      
      f1 <- pheno[idx[1:myf], ]
      f2 <- pheno[idx[(myf+1):(2*myf)], ]
      f3 <- pheno[idx[(2*myf+1):(3*myf)], ]
      f4 <- pheno[idx[(3*myf+1):(4*myf)], ]
      f5 <- pheno[idx[(4*myf+1):length(idx)], ]
      
      
      t1 <- rbind(f2, f3, f4, f5)
      out_t1 <- paste(traits[j], i, "train1.txt", sep="_")
      out_v1 <- paste(traits[j], i, "val1.txt", sep="_")
      write.table(t1, out_t1, sep="\t", row.names=FALSE, quote=FALSE)
      write.table(f1, out_v1, sep="\t", row.names=FALSE, quote=FALSE)
      
      t2 <- rbind(f1, f3, f4, f5)
      out_t2 <- paste(traits[j], i, "train2.txt", sep="_")
      out_v2 <- paste(traits[j], i, "val2.txt", sep="_")
      write.table(t2, out_t2, sep="\t", row.names=FALSE, quote=FALSE)
      write.table(f2, out_v2, sep="\t", row.names=FALSE, quote=FALSE)
      
      t3 <- rbind(f1, f2, f4, f5)
      out_t3 <- paste(traits[j], i, "train3.txt", sep="_")
      out_v3 <- paste(traits[j], i, "val3.txt", sep="_")
      write.table(t3, out_t3, sep="\t", row.names=FALSE, quote=FALSE)
      write.table(f3, out_v3, sep="\t", row.names=FALSE, quote=FALSE)
      
      t4 <- rbind(f1, f2, f3, f5)
      out_t4 <- paste(traits[j], i, "train4.txt", sep="_")
      out_v4 <- paste(traits[j], i, "val4.txt", sep="_")
      write.table(t4, out_t4, sep="\t", row.names=FALSE, quote=FALSE)
      write.table(f4, out_v4, sep="\t", row.names=FALSE, quote=FALSE)
      
      t5 <- rbind(f1, f2, f3, f4)
      out_t5 <- paste(traits[j], i, "train5.txt", sep="_")
      out_v5 <- paste(traits[j], i, "val5.txt", sep="_")
      write.table(t5, out_t5, sep="\t", row.names=FALSE, quote=FALSE)
      write.table(f5, out_v5, sep="\t", row.names=FALSE, quote=FALSE)
      
      message(sprintf("[ %s ] 5 fold split", traits[j]))
    }
  }
  
}  





