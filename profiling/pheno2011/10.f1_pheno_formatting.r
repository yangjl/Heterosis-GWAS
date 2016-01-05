# Jinliang Yang
# Purpose: Formatting BxNAM phenotypic data
# date: 3.26.2011
# location: server.7

setwd("~/Documents/Heterosis_GWAS/pheno2011")
########################################################################################
###### Selected NAM RIL

krn <- read.csv("nam_F1_KRN.csv");
kernel <- read.csv("nam_F1_Kernel.csv")
cob <- read.csv("nam_F1_cob.csv")

trait <- merge(krn, cob, by="Genotype", all=TRUE);
trait <- merge(trait, kernel, by="Genotype", all=TRUE)

trait <- trait[!duplicated(trait$Genotype),]

namril <- subset(trait, Note=="RIL")
namril$Fix <- 1
namril$FID <- NA;
namril$RIL <- as.character(namril$RIL)
for(i in 1:nrow(namril)){
  tem <- unlist(strsplit(namril$RIL[i], split="E"))
  namril$FID[i] <- tem[1]
}
genotypes2 <- namril$RIL
length(genotypes2) #449;

Bxnamril <- subset(trait, Note=="Bx")
Bxnamril$Genotype <- gsub("B73 x ", "", Bxnamril$Genotype)
Bxnamril$Genotype <- gsub(" x B73", "", Bxnamril$Genotype)
Bxnamril$Fix <- 1
Bxnamril$FID <- NA;
for(i in 1:nrow(Bxnamril)){
  tem <- unlist(strsplit(Bxnamril$Genotype[i], split="E"))
  Bxnamril$FID[i] <- tem[1]
}
Bxnamril <- subset(Bxnamril, Genotype%in%genotypes2)
dim(Bxnamril)
#445 16

##########################
cw <- namril[,c("RIL","CW", "Fix", "FID")]
write.table(cw, "pheno_NAMRIL_CW_032612.txt", sep="\t", row.names=FALSE, quote=FALSE)

bxcw <- Bxnamril[,c("RIL","CW", "Fix", "FID")]
write.table(bxcw, "pheno_BxNAMRIL_CW_032612.txt", sep="\t", row.names=FALSE, quote=FALSE)

##########################
cd <- namril[,c("RIL","CD", "Fix", "FID")]
write.table(cd, "pheno_NAMRIL_CD_032612.txt", sep="\t", row.names=FALSE, quote=FALSE)

bxcd <- Bxnamril[,c("RIL","CD", "Fix", "FID")]
write.table(bxcd, "pheno_BxNAMRIL_CD_032612.txt", sep="\t", row.names=FALSE, quote=FALSE)

##########################
cl <- namril[,c("RIL","CL", "Fix", "FID")]
write.table(cl, "pheno_NAMRIL_CL_032612.txt", sep="\t", row.names=FALSE, quote=FALSE)

bxcl <- Bxnamril[,c("RIL","CL", "Fix", "FID")]
write.table(bxcl, "pheno_BxNAMRIL_CL_032612.txt", sep="\t", row.names=FALSE, quote=FALSE)

##########################
kc <- namril[,c("RIL","KC", "Fix", "FID")]
write.table(kc, "pheno_NAMRIL_KC_032612.txt", sep="\t", row.names=FALSE, quote=FALSE)

bxkc <- Bxnamril[,c("RIL","KC", "Fix", "FID")]
write.table(bxkc, "pheno_BxNAMRIL_KC_032612.txt", sep="\t", row.names=FALSE, quote=FALSE)

##########################
akw <- namril[,c("RIL","AKW", "Fix", "FID")]
write.table(akw, "pheno_NAMRIL_AKW_032612.txt", sep="\t", row.names=FALSE, quote=FALSE)

bxakw <- Bxnamril[,c("RIL","AKW", "Fix", "FID")]
write.table(bxakw, "pheno_BxNAMRIL_AKW_032612.txt", sep="\t", row.names=FALSE, quote=FALSE)

##########################
tkw <- namril[,c("RIL","TKW", "Fix", "FID")]
write.table(tkw, "pheno_NAMRIL_TKW_032612.txt", sep="\t", row.names=FALSE, quote=FALSE)

bxtkw <- Bxnamril[,c("RIL","TKW", "Fix", "FID")]
write.table(bxtkw, "pheno_BxNAMRIL_TKW_032612.txt", sep="\t", row.names=FALSE, quote=FALSE)




########################################################################################
###### NAM F1

Bxnamril <- 



Bxkrn <- Bxnamril[, 1:2]


write.table(Bxkrn, "pheno_BxNAMRIL_KRN_030112.txt", sep="\t", row.names=FALSE, quote=FALSE)

Bxkrn2 <- subset(Bxkrn, Genotype %in% namril$RIL)

write.table(Bxkrn2, "pheno_BxNAM445_KRN_032212.txt", sep="\t", row.names=FALSE, quote=FALSE)


genotypes <- Bxnamril$Genotype
genotypes <- gsub("B73 x ", "", genotypes)
genotypes <- gsub(" x B73", "", genotypes)
length(genotypes) #538;

################ tSNP
tSNP <- read.table("~/Documents/Heterosis_GWAS/SNP/tsnp_4impute_AGPv2.txt", header=TRUE)
idx <- which(names(tSNP) %in% genotypes);
dim(tSNP)
#[1] 1055 4897
mytsnp <- tSNP[, c(1:5, idx)]
dim(mytsnp)
#[1] 1055  494
write.table(mytsnp, "tsnp_4BxNAM_AGPv2.txt", sep="\t", row.names=FALSE, quote=FALSE)

##################
./myimpute_beta3.pl --dsnp="fdsnp_2M_012612.dsnp" --tsnp="tsnp_4BxNAM_AGPv2.txt" --start=5 --end=493 --cross="F1" > BxNAM_RIL_GenSel_030112

########################################################################################
###### Selected NAM RIL for training

### find the genotypes for imputation
wmean <- read.table("pheno_KRN_wmean_040211.txt", header=TRUE)
idx2 <- which(wmean$Genotype %in% genotypes);
length(idx2) #527

wmean <- wmean[-idx2, ]
write.table(wmean, "pheno_KRN_wmean_exclude527.txt", quote=FALSE, row.names=FALSE, sep="\t")


################## Impute for PLINK
./myimpute_beta3.pl --dsnp="fdsnp_2M_012612.dsnp" --tsnp="tsnp_4BxNAM_AGPv2.txt" --start=5 --end=493 --cross="F1" --mode=2 --header=2 > BxNAM_RIL_PLINK_030112







