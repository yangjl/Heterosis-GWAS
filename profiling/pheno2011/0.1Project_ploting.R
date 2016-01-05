# Jinliang Yang
# updated: 7/22/2012
# main script for the manipulation and data plotting


# Start to load the project:
setwd("~/Documents/Heterosis_GWAS/pheno2011")
library('ProjectTemplate')
load.project()

############# H2 ###############################
fit1 <- aov(KRN ~ Farm + Genotype, data=nam2011krn)
#16699.7/(16699.7+12575.2+17.1)=0.57 pval=2e-16 ***

fit2 <- aov(CL ~ Farm + Genotype, data=nam2011cob)
#4193402/(4193402+2073463+4091)=0.67 pval=2e-16 ***
fit3 <- aov(CD ~ Farm + Genotype, data=nam2011cob)
#36086/(36086+29667+1233)=0.54 pval=2e-16 ***
fit4 <- aov(CW ~ Farm + Genotype, data=nam2011cob)
#388275/(388275+188621+952)=0.67 pval=2e-16 ***

fit5 <- aov(KC ~ Farm + Genotype, data=nam2011kernel)
#128219105/(128219105+50348796.2+1508212)=0.71 pval=2e-16 ***
fit6 <- aov(TKW ~ Farm + Genotype, data=nam2011kernel)
#12983152/(12983152+4015266+68623)=0.76 pval=2e-16 ***
fit7 <- aov(AKW ~ Farm + Genotype, data=nam2011kernel)
#6.9941/(6.9941+4.2253+0.0149)=0.62 pval=2e-16 ***



#########################################################################################
# plot the HPH
par(mfrow=c(2,4))
hist(diallel$KRN_HPH, main="KRN HPH", xlab="row", cex.lab=1.5);
hist(diallel$KC_HPH, main="KC HPH", xlab="count", cex.lab=1.5);
hist(diallel$TKW_HPH, main="TKW HPH", xlab="weight (g)", cex.lab=1.5);
hist(diallel$AKW_HPH, main="AKW HPH", xlab="weight (g)", cex.lab=1.5);
hist(diallel$CL_HPH, main="CL HPH", xlab="length (mm)", cex.lab=1.5);
hist(diallel$CW_HPH, main="CW HPH", xlab="weight (g)", cex.lab=1.5);
hist(diallel$CD_HPH, main="CD HPH", xlab="diameter (mm)", cex.lab=1.5);

par(mfrow=c(2,4))
hist(diallel$KRN_MPH, main="KRN HPH", xlab="row", cex.lab=2);
hist(diallel$KC_MPH, main="KC HPH", xlab="count", cex.lab=2);
hist(diallel$TKW_MPH, main="TKW HPH", xlab="weight (g)", cex.lab=2);
hist(diallel$AKW_MPH, main="AKW HPH", xlab="weight (g)", cex.lab=2);
hist(diallel$CL_MPH, main="CL HPH", xlab="length (mm)", cex.lab=2);
hist(diallel$CW_MPH, main="CW HPH", xlab="weight (g)", cex.lab=2);
hist(diallel$CD_MPH, main="CD HPH", xlab="diameter (mm)", cex.lab=2);









