# Formatting the 2011 Diallel phenotype data
# Jinliang Yang
# Jan 17th, 2012
# server .7

#setwd("/Users/yangjl/Documents/Heterosis_GWAS/pheno2011")

#diallel <- read.csv("pheno_diallel_master_raw_040612.csv", header=TRUE)
sum(!is.na(diallel$KRN))
# 3225

############# H2 ###############################
fit1 <- aov(KRN ~ Farm + Genotype, data=diallel)
#7974.099/(7974.099+6225.719+221.082)=0.55
fit2 <- aov(KC ~ Farm + Genotype, data=diallel)
#91481212/(91481212+38953037+11996774)=0.64
fit3 <- aov(AKW ~ Farm + Genotype, data=diallel)
#9.2250/(9.2250+4.3802+1.0351)=0.63
fit4 <- aov(TKW ~ Farm + Genotype, data=diallel)
#12480224/(12480224+4092526+1789453)=0.67
fit5 <- aov(CD ~ Farm + Genotype, data=diallel)
#28020.7/(28020.7+17740.9+691.4)=0.60
fit6 <- aov(CW ~ Farm + Genotype, data=diallel)
#407395/(407395+116368+44141)=0.72
fit7 <- aov(CL ~ Farm + Genotype, data=diallel)
#3214417/(1363090+3214417+579579)=0.62


#######################################################
### BOX-COX transformation not necessary
### get the LSmean

library(MASS)
library(ggplot2)

#bct computes the Box-Cox family of transforms: y = (y^lambda - 1)/(lambda*gm^(lambda-1)), 
#where gm is the geometric mean of the y's. returns log(y)*gm when lambda equals 0.
bct <- function(y,lambda){
  gm <- exp( mean( log(y) ) )
  if(lambda==0) return( log(y)*gm )
  yt <- (y^lambda - 1)/( lambda * gm^(lambda-1) )
  return(yt)
}

boxcox(diallel$KRN ~ 1, lambda = seq(-10, 10, 1/10), plotit = TRUE)
boxcox(diallel$KC ~ 1, lambda = seq(-10, 10, 1/10), plotit = TRUE)
boxcox(diallel$AKW ~ 1, lambda = seq(-10, 10, 1/10), plotit = TRUE)
boxcox(diallel$TKW ~ 1, lambda = seq(-10, 10, 1/10), plotit = TRUE)
boxcox(diallel$CD ~ 1, lambda = seq(-10, 10, 1/10), plotit = TRUE)
boxcox(diallel$CW ~ 1, lambda = seq(-10, 10, 1/10), plotit = TRUE)
boxcox(diallel$CL ~ 1, lambda = seq(-10, 10, 1/10), plotit = TRUE)
#ykrn <- bct(diallel$KRN[!is.na(diallel$KRN)], 0)
##################################################################################
library(nlme)
BLUE <- function(data=krn, model=KRN~Genotype2, random=~1|Farm, trait="KRN"){
  lmeout1 <- lme(model, data=data, random=random);
  ped.hat1 <- lmeout1$coef$fixed;
  ped.hat1[-1] <- ped.hat1[-1]+ped.hat1[1];
  names(ped.hat1)[1]="B73xZ001";
  names(ped.hat1) <- gsub("Genotype2", "", names(ped.hat1));
  tped <- data.frame(Genotype=names(ped.hat1), trait=ped.hat1)
  names(tped)[2] <- trait;
    return(tped)
}

###########################################################

# start to calculate BLUE:
krn <- diallel[!is.na(diallel$KRN) & diallel$Type=="Diallel", ]
krn$Note_KRN <- as.character(krn$Note_KRN)
table(krn$Note_KRN)
krn2 <- BLUE(data=krn, model=KRN~Genotype2, trait="KRN")

############### ear trait #####################################
diallel$Note_ear <- as.character(diallel$Note_ear)
table(diallel$Note_ear)
ear <- subset(diallel, KC>50 & Note_ear!=" ")
dim(ear)
# 3057 22
###########-----KC
kc <- ear[!is.na(ear$KC) & ear$Type=="Diallel",]
kc2 <- BLUE(data=kc, model=KC~Genotype2, trait="KC")
hist(kc2$KC)
###########-----AKW
akw <- ear[!is.na(ear$AKW) & ear$Type=="Diallel",]
akw2 <- BLUE(data=akw, model=AKW~Genotype2, trait="AKW")
hist(akw2$AKW)
###########-----TKW
tkw <- ear[!is.na(ear$TKW) & ear$Type=="Diallel",]
tkw2 <- BLUE(data=tkw, model=TKW~Genotype2, trait="TKW")
hist(tkw2$TKW)


############## cob trait #####################################
diallel$Note_cob <- as.character(diallel$Note_cob)
table(diallel$Note_cob)

###########-----CD
cd <- diallel[!is.na(diallel$CD) & diallel$Type=="Diallel",]
cd2 <- BLUE(data=cd, model=CD~Genotype2, trait="CD")
hist(cd2$CD, main="CD")
###########-----CW
cw <- diallel[!is.na(diallel$CW) & diallel$Type=="Diallel",]
cw2 <- BLUE(data=cw, model=CW~Genotype2, trait="CW")
hist(cw2$CW, main="CW")
###########-----CL
cl <- diallel[!is.na(diallel$CL) & diallel$Type=="Diallel",]
cl2 <- BLUE(data=cl, model=CL~Genotype2, trait="CL")
hist(cl2$CL, main="CL")

############################################################################
pheno <- merge(krn2, akw2, by="Genotype", all=TRUE)
pheno <- merge(pheno, kc2, by="Genotype", all=TRUE)
pheno <- merge(pheno, tkw2, by="Genotype", all=TRUE)
pheno <- merge(pheno, cd2, by="Genotype", all=TRUE)
pheno <- merge(pheno, cw2, by="Genotype", all=TRUE)
pheno <- merge(pheno, cl2, by="Genotype", all=TRUE)

write.table(pheno, "cache/pheno_diallel_master_BLUE.csv", sep=",", row.names=FALSE, quote=FALSE)



