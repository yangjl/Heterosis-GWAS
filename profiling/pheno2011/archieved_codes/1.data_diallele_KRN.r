# Formatting the 2011 Diallel phenotype data
# Jinliang Yang
# Jan 17th, 2012
# server .7

setwd("/Users/yangjl/Documents/Heterosis_GWAS/pheno2011")

#############################################################
# Check the data completeness and correct the notes
#Notes: 
#Dairy 8651-9000 v
#Zumwalt 9251-9600 v
#Johnson J2601-J2950 v
#Ki11: kill@ => Ki11 v
#Ms71: MS71 => Ms71 v

file1 <- read.csv("LL_KRN_test.csv")
dim(file1)
#8678 8

diallel <- subset(file1, Note.y=="Diallel" | Note.y=="NAM Filler");
dim(diallel)
#[1] 2850    8
summary(diallel)

###########################################################
### Delete duplicated input and missing data
# Get ride of 4 duplicated records
dup <- which(duplicated(diallel$Barcode))
diallel <- diallel[!duplicated(diallel$Barcode),]

# Get ride of the missing records
### checking for outlyers
diallel$KRN <- as.numeric(as.character(diallel$KRN));
diallel[is.na(diallel$KRN),]
diallel <- diallel[!is.na(diallel$KRN),]
hist(diallel$KRN, breaks=10)
dim(diallel)
#[1] 2837    8

############################################################
# GCA and SCA
diallel$p1 <- NA;
diallel$p2 <- NA;
diallel$reciprocal <- NA;
diallel$Genotype <- as.character(diallel$Genotype);
for(i in 1:nrow(diallel)){
  tem <- unlist(strsplit(diallel$Genotype[i], split=" x "));
  tem <- sort(tem)
  diallel$p1[i] <- tem[1];
  diallel$p2[i] <- tem[2];
  diallel$reciprocal[i] <- paste(tem[1], tem[2], sep="_")
}

fit <- lm(KRN ~ p1 + p2 + reciprocal, data=diallel)

############# H2 ###############################
fit <- aov(KRN ~ Farm + Genotype, data=diallel)
#7458/(7458+5245+0.7)

#######################################################
### get the LSmean
library(nlme)

# compare the difference of the two models and simple mean
lmeout1 <- lme(KRN~Genotype, data=diallel, random=~1|Farm);
ped.hat1 <- lmeout1$coef$fixed;
ped.hat1[-1] <- ped.hat1[-1]+ped.hat1[1];
names(ped.hat1)[1]="B73 x B97";
names(ped.hat1) <- gsub("Genotype","",names(ped.hat1));
tped <- data.frame(Genotype=names(ped.hat1), KRN=ped.hat1)
dim(tped)
#[1] 315   2
###########################################################
### merge the trait of reciprocal F1
library(ggplot2)

tped$Genotype <- as.character(tped$Genotype);
tped$p1 <- NA;
tped$p2 <- NA;
for (i in 1:nrow(tped)){
  tped$p1[i] <- sort(unlist(strsplit(tped$Genotype[i], split=" x ")))[1];
  tped$p2[i] <- sort(unlist(strsplit(tped$Genotype[i], split=" x ")))[2];
}

tped <- tped[order(tped$p1, tped$p2),]
tped$ped <- paste(tped$p1, tped$p2, sep="|")
lsmean_diallel <- ddply(tped, .(ped), summarise,
            KRN=mean(KRN))
dim(lsmean_diallel)
#[1] 248   2

########################################################
# F1 and parents lsmean:
idx <- grep("@", lsmean_diallel$ped)
lsmean <- lsmean_diallel[-idx,]

plsmean <- lsmean_diallel[idx,]
plsmean$ped <- gsub('@\\|NA',"", plsmean$ped)

########################################################
panzeakrn <- read.csv("panzea_KRN.csv")
panzeakrn$Phenotype.Value <- gsub("mean=", "", panzeakrn$Phenotype.Value)

nam <- read.table("nam_parents", header=TRUE)
p <- c("B73", as.character(nam$parent))

panzeakrn <- subset(panzeakrn, Accession %in% p)
panzeakrn$Phenotype.Value <- as.numeric(panzeakrn$Phenotype.Value);
pkrn <- ddply(panzeakrn, .(Accession), summarise,
              KRN=mean(Phenotype.Value, na.rm=T))
pkrn$Accession <- toupper(pkrn$Accession)
plsmean$ped <- toupper(plsmean$ped)
pkrn <- merge(pkrn, plsmean, by.x="Accession", by.y="ped", all=T)

mymean <- function(x) mean(x, na.rm=T);
pkrn$KRN <- apply(pkrn[,2:3], 1, mymean)
pkrn <- pkrn[, c(1,4)]
names(pkrn)[1] <- "ped"

lsmean <- rbind(lsmean, pkrn)
dim(lsmean)
#[1] 252   2
####################################################
## plot and output lsmean
idx2 <- grep("\\|", lsmean$ped);
f1 <- lsmean[idx2,]
p <- lsmean[-idx2, ]
plot(density(f1$KRN), main="KRN of Diallel and parents", xlim=c(8,21), xlab="KRN", lwd=3, col="blue")
lines(density(p$KRN), lwd=3, lty=2, col="red")

write.table(lsmean, "diallel_lsmean_KRN.csv", sep=",", row.names=FALSE, quote=FALSE)

