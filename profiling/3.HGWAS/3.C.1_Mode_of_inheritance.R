# Jinliang Yang
# 8/15/2014
# determining mode of inheritance


#### format genotype
chr1 <- read.csv("~/Documents/Heterosis_GWAS/Method/genoTAV/chr1_recoded.csv", header=FALSE)
chr2 <- read.csv("~/Documents/Heterosis_GWAS/Method/genoTAV/chr2_recoded.csv", header=FALSE)
chr3 <- read.csv("~/Documents/Heterosis_GWAS/Method/genoTAV/chr3_recoded.csv", header=FALSE)
chr4 <- read.csv("~/Documents/Heterosis_GWAS/Method/genoTAV/chr4_recoded.csv", header=FALSE)
chr5 <- read.csv("~/Documents/Heterosis_GWAS/Method/genoTAV/chr5_recoded.csv", header=FALSE)
chr6 <- read.csv("~/Documents/Heterosis_GWAS/Method/genoTAV/chr6_recoded.csv", header=FALSE)
chr7 <- read.csv("~/Documents/Heterosis_GWAS/Method/genoTAV/chr7_recoded.csv", header=FALSE)
chr8 <- read.csv("~/Documents/Heterosis_GWAS/Method/genoTAV/chr8_recoded.csv", header=FALSE)
chr9 <- read.csv("~/Documents/Heterosis_GWAS/Method/genoTAV/chr9_recoded.csv", header=FALSE)
chr10 <- read.csv("~/Documents/Heterosis_GWAS/Method/genoTAV/chr10_recoded.csv", header=FALSE)

#chrall <- rbind(chr5, chr6, chr7, chr8, chr9, chr10)
chrall <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10)
### Note: raw mode, but transformed from PLINK->SNPTEST->raw.
write.table(chrall, "~/Documents/Heterosis_GWAS/HGWAS_proj/cache/geno_chrall_recoded.csv",
            sep=",", row.names=FALSE, quote=FALSE)

######################


ob <- load(file="~/Documents/Heterosis_GWAS/HGWAS_proj/cache/TAVs_7traits.RData")


cdsnpid <- unique(c(lscd[[1]]$snpid, lscd[[2]]$ID, lscd[[3]]$ID))


source("~/Documents/Heterosis_GWAS/HGWAS_proj/lib/getMoI.R")
out1 <- getMoI(ids=lscd[[3]]$ID, trait="CD", fitformula="CD ~ FID + SNP",
               phenofile="~/Documents/Heterosis_GWAS/Method/SNPTEST/CD_chr1-covars.sample")

out2 <- getMoI(ids=lsakw[[3]]$ID, trait="AKW", fitformula="AKW ~ FID + SNP",
               phenofile="~/Documents/Heterosis_GWAS/Method/SNPTEST/AKW_chr1-covars.sample")

out3 <- getMoI(ids=lscl[[3]]$ID, trait="CL", fitformula="CL ~ FID + SNP",
               phenofile="~/Documents/Heterosis_GWAS/Method/SNPTEST/CL_chr1-covars.sample")

out4 <- getMoI(ids=lscw[[3]]$ID, trait="CW", fitformula="CW ~ FID + SNP",
               phenofile="~/Documents/Heterosis_GWAS/Method/SNPTEST/CW_chr1-covars.sample")

out5 <- getMoI(ids=lskc[[3]]$ID, trait="KC", fitformula="KC ~ FID + SNP",
               phenofile="~/Documents/Heterosis_GWAS/Method/SNPTEST/KC_chr1-covars.sample")

out6 <- getMoI(ids=lstkw[[3]]$ID, trait="TKW", fitformula="TKW ~ FID + SNP",
               phenofile="~/Documents/Heterosis_GWAS/Method/SNPTEST/TKW_chr1-covars.sample")


par(mfrow=c(2,3))
plot(c(3,2,1), out1[1, c("AA", "AB", "BB")], xlim=c(1,3), ylim=c(-2, 2), type="l", main="CD")
for(i in 1:nrow(out1)){
  lines(c(3,2,1), out1[i, c("AA", "AB", "BB")])
}
plot(c(3,2,1), out2[1, c("AA", "AB", "BB")], xlim=c(1,3), ylim=c(-2, 2), type="l", main="AKW")
for(i in 1:nrow(out2)){
  lines(c(3,2,1), out2[i, c("AA", "AB", "BB")])
}
plot(c(3,2,1), out3[1, c("AA", "AB", "BB")], xlim=c(1,3), ylim=c(-2, 2), type="l", main="CL")
for(i in 1:nrow(out3)){
  lines(c(3,2,1), out3[i, c("AA", "AB", "BB")])
}
plot(c(3,2,1), out4[1, c("AA", "AB", "BB")], xlim=c(1,3), ylim=c(-2, 2), type="l", main="CW")
for(i in 1:nrow(out4)){
  lines(c(3,2,1), out4[i, c("AA", "AB", "BB")])
}
plot(c(3,2,1), out5[1, c("AA", "AB", "BB")], xlim=c(1,3), ylim=c(-2, 2), type="l", main="KC")
for(i in 1:nrow(out5)){
  lines(c(3,2,1), out5[i, c("AA", "AB", "BB")])
}
plot(c(3,2,1), out6[1, c("AA", "AB", "BB")], xlim=c(1,3), ylim=c(-2, 2), type="l", main="TKW")
for(i in 1:nrow(out6)){
  lines(c(3,2,1), out6[i, c("AA", "AB", "BB")])
}




nrow(subset(out5, h2 < -0.6))

