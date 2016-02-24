# Jinliang Yang
# 5 fold traits split
# date: June/20/2014

### trait split
traits <- c("krn", "cd", "cl", "cw", "akw", "tkw", "kc")
phenopwd <- "/home/NSF-SAM-GS/HGWAS/SNPdiallel/dcv2/pheno"

splitPheno2(traits, splittimes=20, outputpwd=phenopwd, myf=3)


############################  
splitPheno2 <- function(traits, splittimes=10, outputpwd="~", myf=3){
  
  mypheno <- paste("/mnt/02/yangjl/Documents/Heterosis_GWAS/HGWAS_proj/data/", "krn",
                   "_v2_GenSel_diallel.txt", sep="")
  pheno <- read.table(mypheno, header=TRUE)
  pheno$p1 <- gsub("x.*", "", pheno$Genotype)
  pheno$p2 <- gsub(".*x", "", pheno$Genotype)
  
  ps <- unique(c(pheno$p1, pheno$p2))[-1]
  
  setwd(outputpwd)
  for(i in 1:splittimes){
    idx <- sample(1:length(ps), length(ps))
    #myf <- round(length(ps)/7, 0)
    
    for(j in 1:length(traits)){
      mypheno <- paste("/mnt/02/yangjl/Documents/Heterosis_GWAS/HGWAS_proj/data/", traits[j],
                       "_v2_GenSel_diallel.txt", sep="")
      pheno <- read.table(mypheno, header=TRUE)
      pheno$p1 <- gsub("x.*", "", pheno$Genotype)
      pheno$p2 <- gsub(".*x", "", pheno$Genotype)
      
      f1 <- ps[idx[1:myf] ]
      f2 <- ps[idx[(myf+1):(2*myf)]]
      f3 <- ps[idx[(2*myf+1):(3*myf)]]
      f4 <- ps[idx[(3*myf+1):(4*myf)] ]
      f5 <- ps[idx[(4*myf+1):(5*myf)]]
      #f6 <- ps[idx[(5*myf+1):(6*myf)]]
      
      
      #t1 <- subset(pheno, p1 != f1 & p1 != f1)
      v1 <- subset(pheno, p1 %in% f1 | p2 %in% f1)
      t1 <- subset(pheno, !(Genotype %in% v1$Genotype))
      out_t1 <- paste(traits[j], i, "train_s1.txt", sep="_")
      out_v1 <- paste(traits[j], i, "val_s1.txt", sep="_")
      write.table(t1[, 1:3], out_t1, sep="\t", row.names=FALSE, quote=FALSE)
      write.table(v1[, 1:3], out_v1, sep="\t", row.names=FALSE, quote=FALSE)
      
      v2 <- subset(pheno, p1 %in% f2 | p2 %in% f2)
      t2 <- subset(pheno, !(Genotype %in% v2$Genotype))
      out_t2 <- paste(traits[j], i, "train_s2.txt", sep="_")
      out_v2 <- paste(traits[j], i, "val_s2.txt", sep="_")
      write.table(t2[, 1:3], out_t2, sep="\t", row.names=FALSE, quote=FALSE)
      write.table(v2[, 1:3], out_v2, sep="\t", row.names=FALSE, quote=FALSE)
      
      v3 <- subset(pheno, p1 %in% f3 | p2 %in% f3)
      t3 <- subset(pheno, !(Genotype %in% v3$Genotype))
      out_t3 <- paste(traits[j], i, "train_s3.txt", sep="_")
      out_v3 <- paste(traits[j], i, "val_s3.txt", sep="_")
      write.table(t3[, 1:3], out_t3, sep="\t", row.names=FALSE, quote=FALSE)
      write.table(v3[, 1:3], out_v3, sep="\t", row.names=FALSE, quote=FALSE)
      
      v4 <- subset(pheno, p1 %in% f4 | p2 %in% f4)
      t4 <- subset(pheno, !(Genotype %in% v4$Genotype))
      out_t4 <- paste(traits[j], i, "train_s4.txt", sep="_")
      out_v4 <- paste(traits[j], i, "val_s4.txt", sep="_")
      write.table(t4[, 1:3], out_t4, sep="\t", row.names=FALSE, quote=FALSE)
      write.table(v4[, 1:3], out_v4, sep="\t", row.names=FALSE, quote=FALSE)
      
      v5 <- subset(pheno, p1 %in% f5 | p2 %in% f5)
      t5 <- subset(pheno, !(Genotype %in% v5$Genotype))
      out_t5 <- paste(traits[j], i, "train_s5.txt", sep="_")
      out_v5 <- paste(traits[j], i, "val_s5.txt", sep="_")
      write.table(t5[, 1:3], out_t5, sep="\t", row.names=FALSE, quote=FALSE)
      write.table(v5[, 1:3], out_v5, sep="\t", row.names=FALSE, quote=FALSE)
      
      tot <- nrow(pheno)
      message(sprintf("[ %.2f, %.2f, %.2f, %.2f, %.2f ] 5 fold split for [ %s ]", 
                      nrow(v1)/tot, nrow(v2)/tot, nrow(v3)/tot, 
                      nrow(v4)/tot, nrow(v5)/tot, traits[j] ))
    }
  }
  
}  





