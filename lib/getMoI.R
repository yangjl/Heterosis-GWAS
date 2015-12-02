## Jinliang Yang

getMoI <- function(ids=lskc[[3]]$ID, trait="CD", fitformula="CD ~ FID+SNP",
                   phenofile="~/Documents/Heterosis_GWAS/Method/SNPTEST/CD_chr1-covars.sample"){
  #### genotypic data
  mychr <- subset(chrall, V1 %in% ids)
  message(sprintf("# of SNPs: [ %s ]", nrow(mychr)))
  
  #### phenotypic data
  pheno <- read.table(phenofile, header=TRUE)
  pheno <- pheno[-1, 1:6]
  pheno[, trait] <- as.numeric(as.character(pheno[,trait]))
  pheno[, trait] <- scale(pheno[, trait])
  pheno$pop <- paste("P", pheno$pop, sep="")
  pheno$FID <- paste("F", pheno$FID, sep="")
  message(sprintf("# of pheno: [ %s ]", nrow(pheno)))
  
  message(sprintf("# RUNNING for mode of inheritance ...", nrow(pheno)))
  out <- MoI(pheno=pheno, geno=mychr, fitformula=fitformula)
  return(out)
}

##########
MoI <- function(pheno=pheno, geno=mychr, fitformula="CD ~ FID+SNP"){
  
  fitformula <- as.formula(fitformula)
  out <- sapply(1:nrow(geno),
                
                function(i){
                  ph <- pheno;
                  ph$SNP <- t(geno[i, -1])
                  names(ph)[7] <- "SNP"
                  ph <- subset(ph, SNP != "NN")
                  
                  fit <- glm(fitformula, data=ph)
                  fitco <- coef(fit)
                  
                  snpAA <- fitco[1]
                  snpAB <- fitco['SNPAB'] + fitco[1]
                  snpBB <- fitco['SNPBB'] + fitco[1]
                  
                  a <- abs(snpBB-snpAA)/2
                  d <- snpAB - (snpBB+snpAA)/2
                  h1 <- snpAB/((snpAA+snpBB/2))
                  h2 <- d/a
                  
                  res <- as.vector(c(snpAA, snpAB, snpBB, a, d, h1, h2))
                  return(res)
                })
  out <- as.data.frame(t(out))
  names(out) <- c("AA", "AB", "BB", "add", "dom", "h1", "h2")
  row.names(out) <- geno$V1
  return(out)
}
