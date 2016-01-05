# Formatting the 2011 phenotype data: NAM and F1
# Jinliang Yang
# Jan 19th, 2012
# log: change Genotype [. to ""]

setwd("~/Documents/Heterosis_GWAS/pheno2011")

#############################################################
# Check the data completeness and correct the notes
#Notes: Two replications in Curtiss and Johnson
#Curtiss: 4651-6100
#Johnson J1251-J2350

file1 <- read.csv("data/rawdata/LL_cob.test.csv")
dim(file1)
#[1] 8773   10 | [1] 8286   11

nam <- subset(file1, Note.y!="Diallel" & Note.y!="NAM Filler");
dim(nam)
#[1] 5806   10 | [1] 5467   11 | [1] 5828    8
summary(nam)

###########################################################
### Delete duplicated input and missing data
barcode <- nam[duplicated(nam$Barcode),]$Barcode
nam[nam$Barcode%in%barcode,]
length(barcode)
# Get ride of 2 duplicated records
nam <- nam[!duplicated(nam$Barcode),]

##########################################################
### checking for outlyers
nam$CL<- as.numeric(as.character(nam$CL));
nam$CD <- as.numeric(as.character(nam$CD));
nam$CW <- as.numeric(as.character(nam$CW));

hist(nam$CL)
hist(nam$CD)
hist(nam$CW)

head(nam[order(nam$CL),])
tail(nam[order(nam$CL),])

idx <- which.max(nam$CW)
nam <- nam[-idx,]

#extreme long ear: removed
nam <- nam[nam$Barcode!="11-5635-22 OP",]
nam <- nam[nam$Barcode!="11-5622-22 OP",]

nambackup <- nam
dim(nam)
#[1] 5801   10
#######################################################

idx0 <- grep(" x B73", nam$Genotype)
sub1 <- nam[idx0,]
sub1$Genotype <- gsub(" x B73", "", sub1$Genotype)
sub1$Genotype <- paste("B73 x ", sub1$Genotype, sep="")

nam2011cob <- rbind(nam[-idx0,], sub1)

#######################################################
### BLUE

namcl <- BLUE(data=nam2011cob, model=CL~Genotype, random=~1|Farm, trait="CL", intercept="B73 x Z001E0002")
namcl$Genotype <- gsub("@", "", namcl$Genotype)

namcd <- BLUE(data=nam2011cob, model=CD~Genotype, random=~1|Farm, trait="CD", intercept="B73 x Z001E0002")
namcd$Genotype <- gsub("@", "", namcd$Genotype)

namcw <- BLUE(data=nam2011cob, model=CW~Genotype, random=~1|Farm, trait="CW", intercept="B73 x Z001E0002")
namcw$Genotype <- gsub("@", "", namcw$Genotype)
