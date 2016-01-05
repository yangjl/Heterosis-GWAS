# Formatting the 2011 phenotype data: NAM and F1
# Jinliang Yang
# Jan 19th, 2012

setwd("~/Documents/Heterosis_GWAS/pheno2011")

#############################################################
# Check the data completeness and correct the notes
#Notes: Two replications in Curtiss and Johnson
#Curtiss: 4651-6100
#Johnson J1251-J2350

file1 <- read.csv("data/rawdata/LL_KN_AKW_test.csv")
dim(file1)
#[1] 8284   11

nam <- subset(file1, Note.y!="Diallel" & Note.y!="NAM Filler");
dim(nam)
#[1] 5465   11 | [1] 5828    8
summary(nam)

###########################################################
### Delete duplicated input and missing data
barcode <- nam[duplicated(nam$Barcode),]$Barcode
nam[nam$Barcode%in%barcode,]
length(barcode)
# Get ride of 8 duplicated records
nam <- nam[!duplicated(nam$Barcode),]

##########################################################
### checking for outlyers
nam$KC <- as.numeric(as.character(nam$KC));
nam$TKW <- as.numeric(as.character(nam$TKW));
nam$AKW <- as.numeric(as.character(nam$AKW));

#### remove all the records with notes and KC<50
nambackup <- nam
nam <- nam[nam$KC >= 50,]
nam <- nam[nam$Note.x=="",]
hist(nam$KC)

head(nam[order(nam$KC),], 100)
tail(nam[order(nam$KC),])

dim(nam)
#[1] 4737   11
#######################################################

idx0 <- grep(" x B73", nam$Genotype)
sub1 <- nam[idx0,]
sub1$Genotype <- gsub(" x B73", "", sub1$Genotype)
sub1$Genotype <- paste("B73 x ", sub1$Genotype, sep="")

nam2011ear <- rbind(nam[-idx0,], sub1)

#######################################################
### BLUE
namkc <- BLUE(data=nam2011ear, model=KC~Genotype, random=~1|Farm, trait="KC", intercept="B73 x Z001E0002")
namkc$Genotype <- gsub("@", "", namkc$Genotype)

namtkw <- BLUE(data=nam2011ear, model=TKW~Genotype, random=~1|Farm, trait="TKW", intercept="B73 x Z001E0002")
namtkw$Genotype <- gsub("@", "", namtkw$Genotype)

namakw <- BLUE(data=nam2011ear, model=AKW~Genotype, random=~1|Farm, trait="AKW", intercept="B73 x Z001E0002")
namakw$Genotype <- gsub("@", "", namakw$Genotype)

