# purpose: calculate the heterosis
# Jinliang Yang
# 2.1.2012

setwd("/Users/yangjl/Documents/Heterosis_GWAS/pheno2011");
ptrait <- read.csv("diallel_ptrait.csv")

p <- subset(ptrait, Note=="Founder")
f1 <- subset(ptrait, Note=="F1")

f1$KRN_HPH <- NA
f1$KRN_MPH <- NA
f1$CL_HPH <- NA
f1$CL_MPH <- NA
f1$CW_HPH <- NA
f1$CW_MPH <- NA
f1$CD_HPH <- NA
f1$CD_MPH <- NA
f1$AKW_HPH <- NA
f1$AKW_MPH <- NA
f1$TKW_HPH <- NA
f1$TKW_MPH <- NA
f1$KC_HPH <- NA
f1$KC_MPH <- NA

for (i in 1:nrow(f1)){
  p1 <- as.character(f1$p1[i]);
  p2 <- as.character(f1$p2[i]);
  f1$KRN_HPH[i] <- f1$KRN[i] - max(p[p$ped==p1, ]$KRN, p[p$ped==p2, ]$KRN);
  f1$KRN_MPH[i] <- f1$KRN[i] - mean(c(p[p$ped==p1, ]$KRN, p[p$ped==p2, ]$KRN), na.rm=T);
  f1$CL_HPH[i] <- f1$CL[i] - max(p[p$ped==p1, ]$CL, p[p$ped==p2, ]$CL);
  f1$CL_MPH[i] <- f1$CL[i] - mean(c(p[p$ped==p1, ]$CL, p[p$ped==p2, ]$CL), na.rm=T);
  f1$CW_HPH[i] <- f1$CW[i] - max(p[p$ped==p1, ]$CW, p[p$ped==p2, ]$CW);
  f1$CW_MPH[i] <- f1$CW[i] - mean(c(p[p$ped==p1, ]$CW, p[p$ped==p2, ]$CW), na.rm=T);
  f1$CD_HPH[i] <- f1$CD[i] - max(p[p$ped==p1, ]$CD, p[p$ped==p2, ]$CD);
  f1$CD_MPH[i] <- f1$CD[i] - mean(c(p[p$ped==p1, ]$CD, p[p$ped==p2, ]$CD), na.rm=T);
  f1$AKW_HPH[i] <- f1$AKW[i] - max(p[p$ped==p1, ]$AKW, p[p$ped==p2, ]$AKW);
  f1$AKW_MPH[i] <- f1$AKW[i] - mean(c(p[p$ped==p1, ]$AKW, p[p$ped==p2, ]$AKW), na.rm=T);
  f1$TKW_HPH[i] <- f1$TKW[i] - max(p[p$ped==p1, ]$TKW, p[p$ped==p2, ]$TKW);
  f1$TKW_MPH[i] <- f1$TKW[i] - mean(c(p[p$ped==p1, ]$TKW, p[p$ped==p2, ]$TKW), na.rm=T);
  f1$KC_HPH[i] <- f1$KC[i] - max(p[p$ped==p1, ]$KC, p[p$ped==p2, ]$KC);
  f1$KC_MPH[i] <- f1$KC[i] - mean(c(p[p$ped==p1, ]$KC, p[p$ped==p2, ]$KC), na.rm=T);
}
#########################################################################################
par(mfrow=c(2,4))
hist(f1$KRN_HPH, main="KRN HPH", xlab="row", cex.lab=1.5);
hist(f1$KC_HPH, main="KC HPH", xlab="count", cex.lab=1.5);
hist(f1$TKW_HPH, main="TKW HPH", xlab="weight (g)", cex.lab=1.5);
hist(f1$AKW_HPH, main="AKW HPH", xlab="weight (g)", cex.lab=1.5);
hist(f1$CL_HPH, main="CL HPH", xlab="length (mm)", cex.lab=1.5);
hist(f1$CW_HPH, main="CW HPH", xlab="weight (g)", cex.lab=1.5);
hist(f1$CD_HPH, main="CD HPH", xlab="diameter (mm)", cex.lab=1.5);

par(mfrow=c(2,4))
hist(f1$KRN_MPH, main="KRN HPH", xlab="row", cex.lab=2);
hist(f1$KC_MPH, main="KC HPH", xlab="count", cex.lab=2);
hist(f1$TKW_MPH, main="TKW HPH", xlab="weight (g)", cex.lab=2);
hist(f1$AKW_MPH, main="AKW HPH", xlab="weight (g)", cex.lab=2);
hist(f1$CL_MPH, main="CL HPH", xlab="length (mm)", cex.lab=2);
hist(f1$CW_MPH, main="CW HPH", xlab="weight (g)", cex.lab=2);
hist(f1$CD_MPH, main="CD HPH", xlab="diameter (mm)", cex.lab=2);

#################################################################################
diallel <- f1

diallel$p1 <- toupper(diallel$p1)
diallel$p2 <- toupper(diallel$p2)

parents <- read.table("nam_parents", header=TRUE)
parents$parent <- toupper(parents$parent)
B73 <- data.frame(pop="B73", pedigree="B73", parent="B73")
parents <- rbind(B73, parents)

diallel <- merge(diallel, parents[, c(1,3)], by.x="p1", by.y="parent")
diallel <- merge(diallel, parents[, c(1,3)], by.x="p2", by.y="parent")

diallel$FID <- paste(diallel$pop.x, diallel$pop.y, sep="x")
diallel$IID <- 1
diallel <- diallel[,c("FID", "IID", "KRN", "CL", "CW", "CD", "AKW", "TKW", "KC",
                      "KRN_HPH", "CL_HPH", "CW_HPH", "CD_HPH", "AKW_HPH", "TKW_HPH", "KC_HPH",
                      "KRN_MPH", "CL_MPH", "CW_MPH", "CD_MPH", "AKW_MPH", "TKW_MPH", "KC_MPH")];


####default --missing-phenotype -9
diallel[is.na(diallel)] <- -9

write.table(diallel, "diallel_PLINK_040412.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

########################### Specific Combining Ability##########################################






