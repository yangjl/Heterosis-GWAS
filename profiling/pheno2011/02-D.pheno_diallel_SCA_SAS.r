# Jinliang Yang
# purpose: calculate the GCA and SCA
# updated: 3.7.2012
# .7

#################
setwd("~/Documents/Heterosis_GWAS/pheno2011")
rawdata <- read.csv("pheno_diallel_master_raw_040612.csv")
rawdata <- subset(trait, !is.na(pop2))
# 2461 22

data2SAS <- function(input=trait, mytrait="TKW"){
  tkw <- input[, c("Farm", "pop1", "pop2", "Genotype2","Barcode",mytrait)];
  names(tkw) <- c("block","female", "male","cross","Tree","Height");
  outfile <- paste("pheno_diallel_", mytrait, "_4SAS.csv", sep="");
  write.table(tkw, outfile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
}

#############################
data2SAS(input=rawdata, mytrait="KRN")
data2SAS(input=rawdata, mytrait="KC")
data2SAS(input=rawdata, mytrait="AKW")
data2SAS(input=rawdata, mytrait="TKW")
data2SAS(input=rawdata, mytrait="CD")
data2SAS(input=rawdata, mytrait="CL")
data2SAS(input=rawdata, mytrait="CW")


############################
setwd("/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/temoutput")
sca1 <- read.csv("SCA_AKW_BLUP.csv")
sca2 <- read.csv("SCA_KC_BLUP.csv")
sca3 <- read.csv("SCA_KRN_BLUP.csv")
sca4 <- read.csv("SCA_TKW_BLUP.csv")
sca5 <- read.csv("SCA_CD_BLUP.csv")
sca6 <- read.csv("SCA_CW_BLUP.csv")
sca7 <- read.csv("SCA_CL_BLUP.csv")

scatot <- merge(sca1[,2:3], sca2[, 2:3], by="cross")
scatot <- merge(scatot, sca3[, 2:3], by="cross")
scatot <- merge(scatot, sca4[, 2:3], by="cross")
scatot <- merge(scatot, sca5[, 2:3], by="cross")
scatot <- merge(scatot, sca6[, 2:3], by="cross")
scatot <- merge(scatot, sca7[, 2:3], by="cross")

idx <- grep("x", scatot$cross)

SCA <- scatot[idx,]
GCA <- scatot[-idx,]

names(GCA) <- c("Genotype", "AKW", "KC", "KRN", "TKW", "CD", "CW", "CL")
names(SCA) <- c("Genotype", "AKW", "KC", "KRN", "TKW", "CD", "CW", "CL")

setwd("/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/")
write.table(GCA, "cache/pheno_diallel_master_BLUP_GCA.csv", sep=",", row.names=FALSE, quote=FALSE)
write.table(SCA, "cache/pheno_diallel_master_BLUP_SCA.csv", sep=",", row.names=FALSE, quote=FALSE)

