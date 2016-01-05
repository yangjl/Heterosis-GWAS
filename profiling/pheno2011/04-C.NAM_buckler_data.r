# Jinliang Yang
# 7/24/2012
# 


###############################################################################
# NAM data2 from Buckler and Panzea, CW, AKW TKW and KC
###############################################################################

pheno <- read.table("/Users/yangjl/Documents/GWAS2_KRN/pheno/traitMatrix_maize282NAM_v10-120314.txt", header=T)

### Keep only the following traits
### EAR:AKW, KC, TKW=earWeight -CobWeight
### COB: CL

mynm <- names(pheno)
mynm1 <- c("Trait","X20KernelWeight","CobDiameter","CobDiameter.1","CobWeight","CobWeight.1","EarDiameter","EarDiameter.1","EarLength","EarLength.1",
           "EarRowNumber","EarRowNumber.1","EarWeight","EarWeight.1","TotalKernelVolume", "TotalKernelVolume.1")
mypheno <- pheno[, mynm1]
nampheno <- mypheno[-1:-283,]
nampheno[nampheno==-999] <- NA;
nampheno$FID <- NA;
names(nampheno)[1] <- "Genotype";
nampheno$Genotype <- as.character(nampheno$Genotype);
for(i in 1:nrow(nampheno)){
  tem <- unlist(strsplit(nampheno$Genotype[i], split="E"));
  nampheno$FID[i] <- tem[1];
}
nampheno$Fix <- 1;

nampheno[,2:16] <- apply(nampheno[,2:16], 2, as.numeric)
summary(nampheno)

##### AKW #######################################################
nampheno$X20KernelWeight <- as.numeric(as.character(nampheno$X20KernelWeight))/20
names(nampheno)[2] <- "AKW"

##### KRN =EarRowNumber EarRowNumber.1 ###########################


##### KC =TotalKernelVolume TotalKernelVolume.1#######################################################

##### TKW #######################################################
# TKW= EarWeight - CobWeight
nampheno$TKW_rep1 <- nampheno$EarWeight - nampheno$CobWeight
nampheno$TKW_rep2 <- nampheno$EarWeight.1 - nampheno$CobWeight.1

names(nampheno) <- c("Genotype", "AKW","CD_rep1", "CD_rep2", "CW_rep1", "CW_rep2", "EarDiameter",        
 					 "EarDiameter.1", "CL_rep1","CL_rep2", "KRN_rep1","KRN_rep2", "EarWeight", "EarWeight.1",       
					 "KC_rep1", "KC_rep2", "FID","Fix", "TKW_rep1", "TKW_rep2")
library(reshape)
bnamlong <- melt(nampheno, id=c("Genotype", "FID", "Fix"))


### get all the BLUEs
##### CW #######################################################
bcw <- subset(bnamlong, variable=="CW_rep1" | variable=="CW_rep2")
bcw <- bcw[, c(1,4,5)]
names(bcw) <- c("Genotype", "Location", "CW")
bcw <- bcw[bcw$CW<50,]

idx1 <- grep("B73", nam2011blue$Genotype)
nam2011cw <- nam2011blue[-idx1,]
nam2011cw <- nam2011cw[, c("Genotype", "CW")]
nam2011cw$Location <- "Ames2011"
nam2011cw <- nam2011cw[, c(1,3,2)]

cw2008 <- pheno2008[['cw']]
cw2008 <- cw2008[, c("Genotype", "Location", "RIL")]
names(cw2008) <- c("Genotype", "Location", "CW")

hist(bcw$CW)
hist(nam2011cw$CW)
hist(cw2008$CW)

cw <- rbind(bcw, nam2011cw, cw2008)
cw <- cw[!is.na(cw$CW),]
cw_blue <- BLUE(data=cw, model=CW~Genotype, random=~1|Location, trait="CW", intercept="Z001E0001")

##### AKW #######################################################
bakw <- subset(bnamlong, variable=="AKW")
bakw <- bakw[, c(1,4,5)]
names(bakw) <- c("Genotype", "Location", "AKW")

idx1 <- grep("B73", nam2011blue$Genotype)
nam2011akw <- nam2011blue[-idx1,]
nam2011akw <- nam2011akw[, c("Genotype", "AKW")]
nam2011akw$Location <- "Ames2011"
nam2011akw <- nam2011akw[, c(1,3,2)]

akw2008 <- pheno2008[['akw']]
akw2008 <- akw2008[, c("Genotype", "Location", "RIL")]
names(akw2008) <- c("Genotype", "Location", "AKW")

akw <- rbind(bakw, nam2011akw, akw2008)
akw <- akw[!is.na(akw$AKW),]
akw_blue <- BLUE(data=akw, model=AKW~Genotype, random=~1|Location, trait="AKW", intercept="Z001E0001")

##### TKW #######################################################
btkw <- subset(bnamlong, variable=="TKW_rep1" | variable=="TKW_rep2")
btkw <- btkw[, c(1,4,5)]
names(btkw) <- c("Genotype", "Location", "TKW")
btkw <- btkw[btkw$TKW>0,]

idx1 <- grep("B73", nam2011blue$Genotype)
nam2011tkw <- nam2011blue[-idx1,]
nam2011tkw <- nam2011tkw[, c("Genotype", "TKW")]
nam2011tkw$Location <- "Ames2011"
nam2011tkw <- nam2011tkw[, c(1,3,2)]

tkw2008 <- pheno2008[['tkw']]
tkw2008 <- tkw2008[, c("Genotype", "Location", "RIL")]
names(tkw2008) <- c("Genotype", "Location", "TKW")

tkw <- rbind(btkw, nam2011tkw, tkw2008)
tkw <- tkw[!is.na(tkw$TKW),]
tkw_blue <- BLUE(data=tkw, model=TKW~Genotype, random=~1|Location, trait="TKW", intercept="Z001E0001")

##### KC #######################################################
bkc <- subset(bnamlong, variable=="KC_rep1" | variable=="KC_rep2")
bkc <- bkc[, c(1,4,5)]
names(bkc) <- c("Genotype", "Location", "KC")

idx1 <- grep("B73", nam2011blue$Genotype)
nam2011kc <- nam2011blue[-idx1,]
nam2011kc <- nam2011kc[, c("Genotype", "KC")]
nam2011kc$Location <- "Ames2011"
nam2011kc <- nam2011kc[, c(1,3,2)]

kc2008 <- pheno2008[['kc']]
kc2008 <- kc2008[, c("Genotype", "Location", "RIL")]
names(kc2008) <- c("Genotype", "Location", "KC")

kc <- rbind(bkc, nam2011kc, kc2008)
kc <- kc[!is.na(kc$KC),]
kc <- kc[kc$KC>0,]
kc_blue <- BLUE(data=kc, model=KC~Genotype, random=~1|Location, trait="KC", intercept="Z001E0001")
kc_blue <- kc_blue[kc_blue$KC>30,]

#########################
nam_blue <- list()
nam_blue[['akw']] <- akw_blue;
nam_blue[['tkw']] <- tkw_blue;
nam_blue[['kc']] <- kc_blue;
nam_blue[['cw']] <- cw_blue;






