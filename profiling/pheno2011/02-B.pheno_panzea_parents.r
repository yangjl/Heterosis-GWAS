# Formatting the data from panzea phenotype database
# Jinliang Yang
# Jan 17th, 2012
# server .7

setwd("/Users/yangjl/Documents/VariationDB/pheno/panzea_pheno")
library(ggplot2)
############### BLUE function #################################
library(nlme)
BLUE <- function(data=krn, model=KRN~Genotype, random=~1|Datasource/Farm/Year, trait="KRN", intercept="B73"){
  lmeout1 <- lme(model, data=data, random=random);
  ped.hat1 <- lmeout1$coef$fixed;
  ped.hat1[-1] <- ped.hat1[-1]+ped.hat1[1];
  names(ped.hat1)[1]=intercept;
  names(ped.hat1) <- gsub("Genotype", "", names(ped.hat1));
  tped <- data.frame(Genotype=names(ped.hat1), trait=ped.hat1)
  names(tped)[2] <- trait;
  return(tped)
}
######################################################################################
namid <- read.table("~/Documents/Heterosis_GWAS/pheno2011/nam_parents", header=TRUE)
namid <- toupper(namid$parent)

diallel <- read.csv("~/Documents/Heterosis_GWAS/pheno2011/pheno_diallel_master_raw_040612.csv")
idx <- grep("@", diallel$Genotype)
pames <- diallel[idx,]

###### FUNCTION  ##############################################
getblue <- function(amesdata=pameskrn, panzeadata=pkrn, trait="KRN",
                    model=KRN~Genotype, nathan=FALSE){
  ####################### AMES data ----------------------------------
  amesdata <- amesdata[, c("p1", trait, "Farm")];
  names(amesdata)[1] <- "Genotype";
  amesdata[, trait] <- as.numeric(as.character(amesdata[, trait]));
  amesdata <- amesdata[!is.na(amesdata[, trait]),]
    
  blue1 <- BLUE(data=amesdata, model=model, random=~1|Farm, trait=trait);
  
  amesdata$Datasource <- "Ames";
  amesdata$Year <- "2011_Summer"
  
  ####################### Panzea data --------------------------------
  
  if(nathan){
    intercept="B97"
    cols <- c("Env", "INBRED", paste(trait, "Inbred", sep="_"));
    panzeadata <- panzeadata[, cols]
    panzeadata$Datasource <- "Panzea";
    panzeadata$Year <- "2011_Summer";
    panzeadata <- panzeadata[, c("INBRED",paste(trait, "Inbred", sep="_"),"Env","Datasource", "Year")];
    panzeadata <- panzeadata[panzeadata[,paste(trait, "Inbred", sep="_") ]!=".",]
  } else {
    intercept="B73"
    cols <- c("Evaluation.Locality", "Rep", "Accession", "Phenotype.Value");
    panzeadata <- panzeadata[, cols]
    panzeadata$Datasource <- "Panzea";
    panzeadata <- panzeadata[, c("Accession","Phenotype.Value", "Evaluation.Locality", "Datasource", "Rep")];
    panzeadata$Phenotype.Value <- as.character(panzeadata$Phenotype.Value)
    for(i in 1:nrow(panzeadata)){
      tem0 <- unlist(strsplit(panzeadata[i,"Phenotype.Value"], split=";"))
      panzeadata[i,"Phenotype.Value"] <- as.numeric(gsub("mean=", "", tem0[1]))
    }
  }
  
  names(panzeadata) <- names(amesdata)
  panzeadata$Genotype <- toupper(panzeadata$Genotype)
  panzeadata <- subset(panzeadata, Genotype %in% namid)
  panzeadata[, trait] <- as.numeric(as.character(panzeadata[, trait]))
  
  blue2 <- BLUE(data=panzeadata, model=model, random=~1|Datasource/Farm/Year, trait=trait, intercept=intercept)
  
  alldata <- rbind(amesdata, panzeadata)  
  blue3 <- BLUE(data=alldata, model=model, random=~1|Datasource/Farm/Year, trait=trait)
    
  myblue <- merge(blue1, blue2, by="Genotype", all=TRUE)
  myblue<- merge(myblue, blue3, by="Genotype", all=TRUE)
  names(myblue) <- c("Genotype", "Amesdata", "Panzeadata", trait)
  return(myblue)
}

#####################################################################################################
# KRN: ames data, remvoe the annotated ones (normally not well developed);
# KRN: panzea data, (Farm, Year and different datasource)
pames$Note_KRN <- as.character(pames$Note_KRN)
pameskrn <- pames[pames$Note_KRN=="", ]
pkrn <- read.csv("panzea_KRN.csv")
#nathan <- read.csv("~/Documents/Heterosis_GWAS/pheno2011/Nathan_p1_subset.csv")

krnblue <- getblue(amesdata=pameskrn, panzeadata=pkrn, trait="KRN",model=KRN~Genotype, nathan=FALSE)

###### 2. EAR-KC Nathan ##############################################
# KC: ames data, remove the annotated ones and the one have less than 40 observations
# KC: panzea data from Nathan, which only have four locations
pames$Note_ear <- as.character(pames$Note_ear)
ear <- pames[pames$Note_ear=="" & pames$KC >40, ]
nathan <- read.csv("~/Documents/Heterosis_GWAS/pheno2011/Nathan_p1_subset.csv")

kcblue <- getblue(amesdata=ear, panzeadata=nathan, trait="KC", model=KC~Genotype, nathan=TRUE)

###### 3. EAR-TKW Nathan ##############################################
# TKW: ames data, remove the annotated ones and the one have less than 40 observations
# TKW: panzea data from Nathan, which only have four locations
pames$Note_ear <- as.character(pames$Note_ear)
ear <- pames[pames$Note_ear=="" & pames$KC >40, ]
nathan <- read.csv("~/Documents/Heterosis_GWAS/pheno2011/Nathan_p1_subset.csv")

tkwblue <- getblue(amesdata=ear, panzeadata=nathan, trait="TKW", model=TKW~Genotype, nathan=TRUE)

###### 4. EAR-AKW Nathan ##############################################
# TKW: ames data, remove the annotated ones and the one have less than 40 observations
# TKW: panzea data from Nathan, which only have four locations
pames$Note_ear <- as.character(pames$Note_ear)
ear <- pames[pames$Note_ear=="" & pames$KC >40, ]
nathan <- read.csv("~/Documents/Heterosis_GWAS/pheno2011/Nathan_p1_subset.csv")
nathan <- nathan[nathan$X10KW_Inbred !=".",]
nathan$X10KW_Inbred <- as.numeric(as.character(nathan$X10KW_Inbred))/10
names(nathan)[7] <- "AKW_Inbred" 

akwblue <- getblue(amesdata=ear, panzeadata=nathan, trait="AKW", model=AKW~Genotype, nathan=TRUE)

###### 5. COB-CD Nathan ##############################################
# CD: ames data, remove the annotated ones
# CD: panzea data from Nathan, which only have four locations
pames$Note_cob <- as.character(pames$Note_cob)
cob <- pames[pames$Note_cob=="", ]
#nathan <- read.csv("~/Documents/Heterosis_GWAS/pheno2011/Nathan_p1_subset.csv")
pcd <- read.csv("panzea_CD.csv")

cdblue <- getblue(amesdata=cob, panzeadata=pcd, trait="CD", model=CD~Genotype, nathan=FALSE)

###### 6. COB-CW Nathan ##############################################
# CW: ames data, remove the annotated ones
# CW: panzea data from Nathan, which only have four locations
pames$Note_cob <- as.character(pames$Note_cob)
cob <- pames[pames$Note_cob=="", ]
#nathan <- read.csv("~/Documents/Heterosis_GWAS/pheno2011/Nathan_p1_subset.csv")
pcw <- read.csv("panzea_CW.csv")
cwblue <- getblue(amesdata=cob, panzeadata=pcw, trait="CW", model=CW~Genotype, nathan=FALSE)

###### 7. COB-CL Nathan ##############################################
# CL: ames data, remove the annotated ones
# CL: panzea data from Nathan, which only have four locations
pames$Note_cob <- as.character(pames$Note_cob)
cob <- pames[pames$Note_cob=="", ]
nathan <- read.csv("~/Documents/Heterosis_GWAS/pheno2011/Nathan_p1_subset.csv")
#pcl <- read.csv("panzea_KRN.csv")
clblue <- getblue(amesdata=cob, panzeadata=nathan, trait="CL", model=CL~Genotype, nathan=TRUE)
#######################################################################################################################

ptrait <- merge(krnblue[, c(1,4)], kcblue[, c(1,4)], by="Genotype")
ptrait <- merge(ptrait, tkwblue[, c(1,4)], by="Genotype")
ptrait <- merge(ptrait, akwblue[, c(1,4)], by="Genotype")
ptrait <- merge(ptrait, cdblue[, c(1,4)], by="Genotype")
ptrait <- merge(ptrait, cwblue[, c(1,4)], by="Genotype")
ptrait <- merge(ptrait, clblue[, c(1,4)], by="Genotype")

write.table(ptrait, "/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/cache/pheno_namparents_BLUE.csv", sep=",", row.names=FALSE, quote=FALSE)
