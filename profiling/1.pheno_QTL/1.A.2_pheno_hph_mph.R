# Jinliang Yang
# July 18th, 2014
# purpose: get HPH and MPH


###################################

getnamp <- function(){
  namp <- read.csv("data/pheno_namparents_BLUE.csv")
  namid <- read.csv("data/NAM_populations.csv")
  namid$parent <- toupper(namid$parent)
  namp <- merge(namid[, c("parent", "zid")], namp, by.x="parent", by.y="Genotype", all.y=TRUE)
  names(namp)[1] <- "Genotype"
  namp$zid[1] <- "B73"
  return(namp)
}


#################################
source("lib/getHOS.R")

HPH_MPH_of_traits <- function(pheno, namp){
  krn1 <- getHOS(pheno=pheno, namp=namp, mytrait="KRN", mypop="Diallel")
  krn2 <- getHOS(pheno=pheno, namp=namp, mytrait="KRN", mypop="BxRIL")
  krn3 <- getHOS(pheno=pheno, namp=namp, mytrait="KRN", mypop="MxRIL")
  krn <- rbind(krn1, krn2, krn3)
  message(sprintf("avg pHPH [ %.3f ] and pMPH [ %.3f ]", mean(krn$pHPH), mean(krn$pMPH)))
  
  cd1 <- getHOS(pheno=pheno, namp=namp, mytrait="CD", mypop="Diallel")
  cd2 <- getHOS(pheno=pheno, namp=namp, mytrait="CD", mypop="BxRIL")
  cd3 <- getHOS(pheno=pheno, namp=namp, mytrait="CD", mypop="MxRIL")
  cd <- rbind(cd1, cd2, cd3)
  message(sprintf("avg pHPH [ %.3f ] and pMPH [ %.3f ]", mean(cd$pHPH), mean(cd$pMPH)))
  
  akw1 <- getHOS(pheno=pheno, namp=namp, mytrait="AKW", mypop="Diallel")
  akw2 <- getHOS(pheno=pheno, namp=namp, mytrait="AKW", mypop="BxRIL")
  akw3 <- getHOS(pheno=pheno, namp=namp, mytrait="AKW", mypop="MxRIL")
  akw <- rbind(akw1, akw2, akw3)
  message(sprintf("avg pHPH [ %.3f ] and pMPH [ %.3f ]", mean(akw$pHPH), mean(akw$pMPH)))
  
  cl1 <- getHOS(pheno=pheno, namp=namp, mytrait="CL", mypop="Diallel")
  cl2 <- getHOS(pheno=pheno, namp=namp, mytrait="CL", mypop="BxRIL")
  cl3 <- getHOS(pheno=pheno, namp=namp, mytrait="CL", mypop="MxRIL")
  cl <- rbind(cl1, cl2, cl3)
  message(sprintf("avg pHPH [ %.3f ] and pMPH [ %.3f ]", mean(cl$pHPH), mean(cl$pMPH)))
  
  cw1 <- getHOS(pheno=pheno, namp=namp, mytrait="CW", mypop="Diallel")
  cw2 <- getHOS(pheno=pheno, namp=namp, mytrait="CW", mypop="BxRIL")
  cw3 <- getHOS(pheno=pheno, namp=namp, mytrait="CW", mypop="MxRIL")
  cw <- rbind(cw1, cw2, cw3)
  message(sprintf("avg pHPH [ %.3f ] and pMPH [ %.3f ]", mean(cw$pHPH), mean(cw$pMPH)))
  
  kc1 <- getHOS(pheno=pheno, namp=namp, mytrait="KC", mypop="Diallel")
  kc2 <- getHOS(pheno=pheno, namp=namp, mytrait="KC", mypop="BxRIL")
  kc3 <- getHOS(pheno=pheno, namp=namp, mytrait="KC", mypop="MxRIL")
  kc <- rbind(kc1, kc2, kc3)
  message(sprintf("avg pHPH [ %.3f ] and pMPH [ %.3f ]", mean(kc$pHPH), mean(kc$pMPH)))
  
  tkw1 <- getHOS(pheno=pheno, namp=namp, mytrait="TKW", mypop="Diallel")
  tkw2 <- getHOS(pheno=pheno, namp=namp, mytrait="TKW", mypop="BxRIL")
  tkw3 <- getHOS(pheno=pheno, namp=namp, mytrait="TKW", mypop="MxRIL")
  tkw <- rbind(tkw1, tkw2, tkw3)
  message(sprintf("avg pHPH [ %.3f ] and pMPH [ %.3f ]", mean(tkw$pHPH), mean(tkw$pMPH)))
  
  ###
  #source("~/Documents/Github/zmSNPtools/Rcodes/save.append.R")
  save(list=c("krn", "cd", "akw", "cl", "cw", "kc", "tkw"), 
              file="cache/heterosis_traits.RData")
}

#####
pheno <- read.table("reports/S1_pheno.txt", header=TRUE)
namp <- getnamp()
HPH_MPH_of_traits(pheno, namp)



