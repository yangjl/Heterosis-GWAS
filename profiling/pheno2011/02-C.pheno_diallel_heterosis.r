# purpose: calculate the heterosis (MPH and HPH)
# Jinliang Yang
# 2.1.2012
# updated: 4.12.2012

setwd("/Users/yangjl/Documents/Heterosis_GWAS/pheno2011");
#### Heterosis ######################################################

diallel <- read.csv("cache/pheno_diallel_master_BLUE.csv")
dim(diallel)
# 225 8
p <- read.csv("cache/pheno_namparents_blue.csv")
dim(p)
# 27 8

namid <- read.table("nam_parents", header=TRUE)
namid <- namid[, c(1,3)]

diallel$z1 <- NA;
diallel$z2 <- NA;
diallel$Genotype <- as.character(diallel$Genotype)
for(i in 1:nrow(diallel)){
  tem <- unlist(strsplit(diallel$Genotype[i], split="x"))
  diallel$z1[i] <- tem[1];
  diallel$z2[i] <- tem[2];
}

diallel <- merge(namid, diallel, by.x="pop", by.y="z1")
diallel <- merge(namid, diallel, by.x="pop", by.y="z2")
diallel <- diallel[, c(5, 3,1, 4,2, 6:12)]
names(diallel)[2:5] <- c("z1", "z2", "p1", "p2")
diallel$p1 <- toupper(diallel$p1);
diallel$p2 <- toupper(diallel$p2);

diallel$KRN_HPH <- NA
diallel$KRN_MPH <- NA
diallel$CL_HPH <- NA
diallel$CL_MPH <- NA
diallel$CW_HPH <- NA
diallel$CW_MPH <- NA
diallel$CD_HPH <- NA
diallel$CD_MPH <- NA
diallel$AKW_HPH <- NA
diallel$AKW_MPH <- NA
diallel$TKW_HPH <- NA
diallel$TKW_MPH <- NA
diallel$KC_HPH <- NA
diallel$KC_MPH <- NA

##########################################################################################################
p$Genotype <- as.character(p$Genotype)
for (i in 1:nrow(diallel)){
  p1 <- as.character(diallel$p1[i]);
  p2 <- as.character(diallel$p2[i]);
  diallel$KRN_HPH[i] <- diallel$KRN[i] - max(p[p$Genotype==p1, ]$KRN, p[p$Genotype==p2, ]$KRN);
  diallel$KRN_MPH[i] <- diallel$KRN[i] - mean(c(p[p$Genotype==p1, ]$KRN, p[p$Genotype==p2, ]$KRN), na.rm=T);
  diallel$CL_HPH[i] <- diallel$CL[i] - max(p[p$Genotype==p1, ]$CL, p[p$Genotype==p2, ]$CL);
  diallel$CL_MPH[i] <- diallel$CL[i] - mean(c(p[p$Genotype==p1, ]$CL, p[p$Genotype==p2, ]$CL), na.rm=T);
  diallel$CW_HPH[i] <- diallel$CW[i] - max(p[p$Genotype==p1, ]$CW, p[p$Genotype==p2, ]$CW);
  diallel$CW_MPH[i] <- diallel$CW[i] - mean(c(p[p$Genotype==p1, ]$CW, p[p$Genotype==p2, ]$CW), na.rm=T);
  diallel$CD_HPH[i] <- diallel$CD[i] - max(p[p$Genotype==p1, ]$CD, p[p$Genotype==p2, ]$CD);
  diallel$CD_MPH[i] <- diallel$CD[i] - mean(c(p[p$Genotype==p1, ]$CD, p[p$Genotype==p2, ]$CD), na.rm=T);
  diallel$AKW_HPH[i] <- diallel$AKW[i] - max(p[p$Genotype==p1, ]$AKW, p[p$Genotype==p2, ]$AKW);
  diallel$AKW_MPH[i] <- diallel$AKW[i] - mean(c(p[p$Genotype==p1, ]$AKW, p[p$Genotype==p2, ]$AKW), na.rm=T);
  diallel$TKW_HPH[i] <- diallel$TKW[i] - max(p[p$Genotype==p1, ]$TKW, p[p$Genotype==p2, ]$TKW);
  diallel$TKW_MPH[i] <- diallel$TKW[i] - mean(c(p[p$Genotype==p1, ]$TKW, p[p$Genotype==p2, ]$TKW), na.rm=T);
  diallel$KC_HPH[i] <- diallel$KC[i] - max(p[p$Genotype==p1, ]$KC, p[p$Genotype==p2, ]$KC);
  diallel$KC_MPH[i] <- diallel$KC[i] - mean(c(p[p$Genotype==p1, ]$KC, p[p$Genotype==p2, ]$KC), na.rm=T);
}

write.table(diallel, "cache/pheno_diallel_master_BLUE_heterosis.csv", sep=",", quote=FALSE, row.names=FALSE)

