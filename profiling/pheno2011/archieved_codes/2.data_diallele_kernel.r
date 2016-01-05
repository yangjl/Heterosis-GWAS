# Formatting the 2011 Diallel phenotype data
# Jinliang Yang
# Jan 17th, 2012
# KC, TKW and AKW

#setwd("/Users/yangjl/Documents/workingSpace/pheno2011/codes")

setwd("/Users/yangjl/Documents/Heterosis_GWAS/pheno2011")

#############################################################
# Check the data completeness and correct the notes
#Notes: 
#Dairy 8651-9000 v
#Zumwalt 9251-9600 v
#Johnson J2601-J2950 v
#Ki11: kill@ => Ki11 v
#Ms71: MS71 => Ms71 v

file1 <- read.csv("LL_KN_AKW_test.csv")
dim(file1)
#[1] 8286   11 ||||| #8678 8

diallel <- subset(file1, Note.y=="Diallel" | Note.y=="NAM Filler"| Note.y=="B73");
dim(diallel)
#[1] 3204   11 | [1] 2850    8
summary(diallel)

###########################################################
### Delete duplicated input and missing data
# Get ride of three duplicated records
barcode <- diallel[duplicated(diallel$Barcode),]$Barcode
diallel[diallel$Barcode %in% barcode,]
diallel <- diallel[!duplicated(diallel$Barcode),]
dim(diallel)
#[1] 3200   11
# Get ride of recodes with note.x nad KC<5
diallel <- subset(diallel, Note.x=="" & KC >= 50)
#diallel[is.na(diallel$TKW),]
#diallel <- diallel[!is.na(diallel$KRN),]
##########################################################
### checking for outlyers
diallel$KC <- as.numeric(as.character(diallel$KC));
hist(diallel$AKW)
hist(diallel$KC, xlim=c(0, 1000))
hist(diallel$TKW)

#######################################################
### get the LSmean
library(nlme)

#fit the model (it takes a while to run in R)
lmeout1 <- lme(KC~Genotype, data=diallel, random=~1|Farm);
lmeout2 <- lme(TKW~Genotype, data=diallel, random=~1|Farm);
lmeout3 <- lme(AKW~Genotype, data=diallel, random=~1|Farm);

#extract the estimates of coefficents, which are the lsmeans 
ped.hat1 <- lmeout1$coef$fixed
ped.hat2 <- lmeout2$coef$fixed
ped.hat3 <- lmeout3$coef$fixed

#reparameterize the estimates to be more interpretable 
ped.hat1[-1] <- ped.hat1[-1]+ped.hat1[1];
ped.hat2[-1] <- ped.hat2[-1]+ped.hat2[1];
ped.hat3[-1] <- ped.hat3[-1]+ped.hat3[1];

#rename the first estimate
head(ped.hat1)
head(diallel[order(diallel$Genotype),], 30);

names(ped.hat1)[1]="B73 x B97";
names(ped.hat2)[1]="B73 x B97";
names(ped.hat3)[1]="B73 x B97";

#get rid of the prefix "Genotype" in the levels names so that they are more readable
names(ped.hat1) <- gsub("Genotype","",names(ped.hat1));
names(ped.hat2) <- gsub("Genotype","",names(ped.hat2));
names(ped.hat3) <- gsub("Genotype","",names(ped.hat3));
#transpose the results
kc <- data.frame(Genotype=names(ped.hat1), KC=ped.hat1);
tkw <- data.frame(Genotype=names(ped.hat2), TKW=ped.hat2);
akw <- data.frame(Genotype=names(ped.hat3), AKW=ped.hat3);

############# H2 ###############################
fit1 <- aov(KC ~ Farm + Genotype, data=diallel)
#71424477/(71424477+33294773+9504759)=62.5 pval=2.2e-16 ***
fit2 <- aov(TKW ~ Farm + Genotype, data=diallel)
#10406226/(3627354+10406226+1473890)=67.1 pval=2.2e-16 ***
fit3 <- aov(AKW ~ Farm + Genotype, data=diallel)
#8.3488/(8.3488+3.8944+0.8990)=63.5 pval=2.2e-16 ***
###########################################################
### merge the trait of reciprocal F1
library(ggplot2)
########-------KC
kc$Genotype <- as.character(kc$Genotype);
kc$p1 <- NA;
kc$p2 <- NA;
for (i in 1:nrow(kc)){
  tem <- sort(unlist(strsplit(kc$Genotype[i], split=" x ")));
  kc$p1[i] <- tem[1];
  kc$p2[i] <- tem[2];
}

kc <- kc[order(kc$p1, kc$p2),]
kc$ped <- paste(kc$p1, kc$p2, sep="|")
lsmean_kc <- ddply(kc, .(ped), summarise,
                   KC=mean(KC))
dim(lsmean_kc)
#[1] 249   2

########-------TKW
tkw$Genotype <- as.character(tkw$Genotype);
tkw$p1 <- NA;
tkw$p2 <- NA;
for (i in 1:nrow(tkw)){
  tem <- sort(unlist(strsplit(tkw$Genotype[i], split=" x ")));
  tkw$p1[i] <- tem[1];
  tkw$p2[i] <- tem[2];
}

tkw <- tkw[order(tkw$p1, tkw$p2),]
tkw$ped <- paste(tkw$p1, tkw$p2, sep="|")
lsmean_tkw <- ddply(tkw, .(ped), summarise,
                   TKW=mean(TKW))
dim(lsmean_tkw)
#[1] 249   2

########-------AKW
akw$Genotype <- as.character(akw$Genotype);
akw$p1 <- NA;
akw$p2 <- NA;
for (i in 1:nrow(akw)){
  tem <- sort(unlist(strsplit(akw$Genotype[i], split=" x ")));
  akw$p1[i] <- tem[1];
  akw$p2[i] <- tem[2];
}

akw <- akw[order(akw$p1, akw$p2),]
akw$ped <- paste(akw$p1, akw$p2, sep="|")
lsmean_akw <- ddply(akw, .(ped), summarise,
                   AKW=mean(AKW))
dim(lsmean_akw)
#[1] 249   2
######################################################################
nathan <- read.csv("Nathan_p1_subset.csv")
nam <- read.table("nam_parents", header=TRUE)
p <- c("B73", as.character(nam$parent))
p <- toupper(p);
nathan$INBRED <- toupper(nathan$INBRED)

pcob <- subset(nathan, INBRED %in% p)
pcob[pcob=="."] <- NA;
pcob$X10KW_Inbred <- as.numeric(as.character(pcob$X10KW_Inbred))
pcob$TKW_Inbred <- as.numeric(as.character(pcob$TKW_Inbred))
pcob$KC_Inbred <- as.numeric(as.character(pcob$KC_Inbred))
pcobmean <- ddply(pcob, .(INBRED), summarise,
                  AKW=mean(X10KW_Inbred, na.rm=T)/10,
                  TKW=mean(TKW_Inbred, na.rm=T),
                  KC=mean(KC_Inbred, na.rm=T)
                  )
b73 <- data.frame(INBRED="B73", AKW=0.2583333 , TKW=62.85, KC=243.2903)
pcobmean <- rbind(b73, pcobmean)
########################################################
#------------------ AKW F1 and parents lsmean:
idx1 <- grep("@", lsmean_akw$ped)
f1akw <- lsmean_akw[-idx1,]

pakw <- lsmean_akw[idx1,]
pakw$ped <- gsub('@\\|NA',"", pakw$ped)

pakw$ped <- toupper(pakw$ped);
pakw2 <- merge(pakw, pcobmean[, c("INBRED", "AKW")], by.x="ped", by.y="INBRED", all=T)

cor(subset(pakw2, !is.na(AKW.x) & !is.na(AKW.y))[,2:3])
#0.50

mymean <- function(x) mean(x, na.rm=T);
pakw2$AKW <- apply(pakw2[,2:3], 1, mymean)
pakw2 <- pakw2[, c(1,4)]
names(pakw2)[1] <- "ped"

master_akw <- rbind(f1akw, pakw2)
dim(master_akw)
#[1] 254   2
#------------------ TKW F1 and parents lsmean:
idx2 <- grep("@", lsmean_tkw$ped)
f1tkw <- lsmean_tkw[-idx2,]

ptkw <- lsmean_tkw[idx2,]
ptkw$ped <- gsub('@\\|NA',"", ptkw$ped)

ptkw$ped <- toupper(ptkw$ped);
ptkw2 <- merge(ptkw, pcobmean[, c("INBRED", "TKW")], by.x="ped", by.y="INBRED", all=T)

cor(subset(ptkw2, !is.na(TKW.x) & !is.na(TKW.y))[,2:3])
#0.45

#mymean <- function(x) mean(x, na.rm=T);
#pcw2$CW <- apply(pcw2[,2:3], 1, mymean)
ptkw2 <- ptkw2[, c(1,3)]
names(ptkw2)[2] <- "TKW"

master_tkw <- rbind(f1tkw, ptkw2)

#------------------ KC F1 and parents lsmean:
idx3 <- grep("@", lsmean_kc$ped)
f1kc <- lsmean_kc[-idx3,]

pkc <- lsmean_kc[idx3,]
pkc$ped <- gsub('@\\|NA',"", pkc$ped)

pkc$ped <- toupper(pkc$ped);
pkc2 <- merge(pkc, pcobmean[, c("INBRED", "KC")], by.x="ped", by.y="INBRED", all=T)

cor(subset(pkc2, !is.na(KC.x) & !is.na(KC.y))[,2:3])
#0.61

#mymean <- function(x) mean(x, na.rm=T);
#pcd2$CD <- apply(pcd2[,2:3], 1, mymean)
pkc2 <- pkc2[, c(1,3)]
names(pkc2)[2] <- "KC"

master_kc <- rbind(f1kc, pkc2)
###########################################
#### Merge the data for output
master_kernel <- merge(master_akw, master_tkw, by="ped", all=T)
master_kernel <- merge(master_kernel, master_kc, by="ped", all=T)

write.table(master_kernel, "diallel_lsmean_kernel_020112.csv", sep=",", row.names=FALSE, quote=FALSE)
## plot and output lsmean
idx2 <- grep("\\|", master_kc$ped);
f1 <- master_kc[idx2,]
p <- master_kc[-idx2, ]
plot(density(f1$KC), main="KC of Diallel and parents", xlab="KRN", lwd=3, col="blue")
lines(density(p$KC), lwd=3, lty=2, col="red")

