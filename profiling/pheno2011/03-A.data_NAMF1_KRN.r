# Formatting the 2011 phenotype data: NAMRILs and Bx NAMRILs
# Jinliang Yang
# Jan 19th, 2012
# location: .7
# log: change [. to ""]

#############################################################
# Check the data completeness and correct the notes
#Notes: Two replications in Curtiss and Johnson
#Curtiss: 4651-6100
#Johnson J1251-J2350

file1 <- read.csv("~/Documents/workingSpace/Heterosis_GWAS/pheno2011/2011rawdata/LL_KRN_test.csv")
dim(file1)
#8678 8

nam <- subset(file1, Note.y!="Diallel" & Note.y!="NAM Filler");
dim(nam)
#[1] 5828    8
summary(nam)

###########################################################
### Delete duplicated input and missing data
barcode <- nam[duplicated(nam$Barcode),]$Barcode
nam[nam$Barcode%in%barcode,]

# Get ride of 13 duplicated records
nam <- nam[!duplicated(nam$Barcode),]

##########################################################
### checking for outlyers
nam$KRN <- as.numeric(as.character(nam$KRN));
hist(nam$KRN, breaks=30)
# Get ride of the missing records
nrow(nam[is.na(nam$KRN),]);
nam <- nam[!is.na(nam$KRN),]
dim(nam)
#[1] 5790    8

idx0 <- grep(" x B73", nam$Genotype)
sub1 <- nam[idx0,]
sub1$Genotype <- gsub(" x B73", "", sub1$Genotype)
sub1$Genotype <- paste("B73 x ", sub1$Genotype, sep="")

nam2011krn <- rbind(nam[-idx0,], sub1)

#######################################################
### get the LSmean
library(nlme)

namkrn <- BLUE(data=nam2011krn, model=KRN~Genotype, random=~1|Farm, trait="KRN", intercept="B73 x Z001E0002")
namkrn$Genotype <- gsub("@", "", namkrn$Genotype)

