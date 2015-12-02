# Jinliang Yang
# July 22nd, 2014
# calculate MPH and HPH
###################################
getHOS <- function(pheno=pheno, namp=namp, mytrait="KRN", mypop="Diallel"){
  #mypop: BxRIL Diallel   MxRIL  NAMRIL
  
  myp <- subset(pheno, trait %in% mytrait & pop %in% mypop)
  if(mypop == "BxRIL"){
    message(sprintf("BxRIL: [ %s ] phenotype of [ %s ] lines", mytrait, nrow(myp)))
    myp$p1 <- gsub("x.*", "", myp$Genotype)
    myp$p2 <- gsub(".*x", "", myp$Genotype)
    myp$p1t <- namp[namp$Genotype=="B73", mytrait]
    
    tempheno <- subset(pheno, trait %in% mytrait)[, c("Genotype", "value")]
    names(tempheno)[2] <- "p2t"
    myp <- merge(myp, tempheno, by.x="p2", by.y="Genotype")
    message(sprintf("[ %s ] lines have parental info!", nrow(myp) ))
    myp <- calHOS(myp=myp)
  }
  ######
  if(mypop == "MxRIL"){
    message(sprintf("MxRIL: [ %s ] phenotype of [ %s ] lines", mytrait, nrow(myp)))
    myp$p1 <- gsub("x.*", "", myp$Genotype)
    myp$p2 <- gsub(".*x", "", myp$Genotype)
    myp$p1t <- namp[namp$Genotype=="MO17", mytrait]
    
    tempheno <- subset(pheno, trait %in% mytrait)[, c("Genotype", "value")]
    names(tempheno)[2] <- "p2t"
    myp <- merge(myp, tempheno, by.x="p2", by.y="Genotype")
    message(sprintf("[ %s ] lines have parental info!", nrow(myp) ))
    myp <- calHOS(myp=myp)
  }
  #####
  if(mypop == "Diallel"){
    message(sprintf("Diallel: [ %s ] phenotype of [ %s ] lines", mytrait, nrow(myp)))
    myp$p1 <- gsub("x.*", "", myp$Genotype)
    myp$p2 <- gsub(".*x", "", myp$Genotype)
    
    tempheno1 <- namp[, c("zid", mytrait)]
    names(tempheno1)[2] <- "p1t"
    myp <- merge(myp, tempheno1, by.x="p1", by.y="zid")
    tempheno2 <- namp[, c("zid", mytrait)]
    names(tempheno2)[2] <- "p2t"
    myp <- merge(myp, tempheno2, by.x="p2", by.y="zid")
    message(sprintf("[ %s ] lines have parental info!", nrow(myp) ))
    myp <- calHOS(myp=myp)
  }
  message(sprintf("Avg. HPH [%.3f], pHPH [%.3f], MPH [%.3f] and pMPH [%.3f]",
                  mean(myp$HPH), mean(myp$pHPH), mean(myp$MPH), mean(myp$pMPH)))
  return(myp)
}

calHOS <- function(myp=myp){
  #myp: p1t (P1), p2t(P2) and value(F1)
  myp$pMPH <- myp$MPH <- myp$pHPH <- myp$HPH <- 0;
  for(i in 1:nrow(myp)){
    myp$HPH[i] <- myp$value[i] - max(myp$p1t[i], myp$p2t[i])
    myp$pHPH[i] <- (myp$value[i] - max(myp$p1t[i], myp$p2t[i]))/max(myp$p1t[i], myp$p2t[i])
    myp$MPH[i] <- myp$value[i] - mean(myp$p1t[i], myp$p2t[i])
    myp$pMPH[i] <- (myp$value[i] - mean(myp$p1t[i], myp$p2t[i]))/mean(myp$p1t[i], myp$p2t[i])
  }
  return(myp) 
}
