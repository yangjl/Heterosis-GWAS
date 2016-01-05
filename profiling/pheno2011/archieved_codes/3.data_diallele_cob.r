# Formatting the 2011 Diallel phenotype data
# Jinliang Yang
# Jan 17th, 2012
# CL, CD and CW: Diallel Cob Trait

#setwd("/Users/yangjl/Documents/workingSpace/pheno2011/2011rawdata")

setwd("/Users/yangjl/Documents/pheno2011")

#############################################################
# Check the data completeness and correct the notes
#Notes: 
#Dairy 8651-9000 v
#Zumwalt 9251-9600 v
#Johnson J2601-J2950 v
#Ki11: kill@ => Ki11 v
#Ms71: MS71 => Ms71 v

file1 <- read.csv("LL_cob.test.csv")
dim(file1)
#[1] 8774   10| kernel [1] 8286   11 |KRN #8678 8

diallel <- subset(file1, Note.y=="Diallel" | Note.y=="NAM Filler" | Note.y=="B73");
dim(diallel)
#[1] 2966   10 | [1] 2819   11 | [1] 2850    8
summary(diallel)

###########################################################
### Delete duplicated input and missing data
# Get ride of 126 duplicated records
barcode <- diallel[duplicated(diallel$Barcode),]$Barcode
diallel[diallel$Barcode %in% barcode,]
diallel <- diallel[!duplicated(diallel$Barcode),]

# Get ride of the missing records
#diallel[is.na(diallel$TKW),]
#diallel <- diallel[!is.na(diallel$KRN),]
##########################################################
### checking for outlyers
diallel$CL <- as.numeric(as.character(diallel$CL));
diallel$CD <- as.numeric(as.character(diallel$CD));
diallel$CW <- as.numeric(as.character(diallel$CW));


hist(diallel$CL);
hist(diallel$CD);
hist(diallel$CW);
idx <- which.max(diallel$CW)
diallel[idx,]
diallel[idx, ]$CW <- 28.99
#######################################################
### get the LSmean
library(nlme)

#fit the model (it takes a while to run in R)
tem <- diallel[!is.na(diallel$CL),]
lmeout1 <- lme(CL~Genotype, data=tem, random=~1|Farm);
lmeout2 <- lme(CD~Genotype, data=diallel, random=~1|Farm);
lmeout3 <- lme(CW~Genotype, data=diallel, random=~1|Farm);

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
cl <- data.frame(Genotype=names(ped.hat1), CL=ped.hat1);
cd <- data.frame(Genotype=names(ped.hat2), CD=ped.hat2);
cw <- data.frame(Genotype=names(ped.hat3), CW=ped.hat3);

############# H2 ###############################
fit1 <- aov(CL ~ Farm + Genotype, data=diallel)
#2954644/(2954644+39598+1261566)=69.4 pval=2.2e-16 ***
fit2 <- aov(CD ~ Farm + Genotype, data=diallel)
#28022/(28022+15250+174)=64.5 pval=2.2e-16 ***
fit3 <- aov(CW ~ Farm + Genotype, data=diallel)
#390226/(390226+108061+2492)=77.9 pval=2.2e-16 ***
###########################################################

### merge the trait of reciprocal F1
library(ggplot2)
########-------CL
cl$Genotype <- as.character(cl$Genotype);
cl$p1 <- NA;
cl$p2 <- NA;
for (i in 1:nrow(cl)){
  tem <- sort(unlist(strsplit(cl$Genotype[i], split=" x ")));
  cl$p1[i] <- tem[1];
  cl$p2[i] <- tem[2];
}

cl <- cl[order(cl$p1, cl$p2),]
cl$ped <- paste(cl$p1, cl$p2, sep="|")
lsmean_cl <- ddply(cl, .(ped), summarise,
                        CL=mean(CL))
dim(lsmean_cl)
#[1] 251   2

########-------CW
cw$Genotype <- as.character(cw$Genotype);
cw$p1 <- NA;
cw$p2 <- NA;
for (i in 1:nrow(cw)){
  tem <- sort(unlist(strsplit(cw$Genotype[i], split=" x ")));
  cw$p1[i] <- tem[1];
  cw$p2[i] <- tem[2];
}

cw <- cw[order(cw$p1, cw$p2),]
cw$ped <- paste(cw$p1, cw$p2, sep="|")
lsmean_cw <- ddply(cw, .(ped), summarise,
                   CW=mean(CW))
dim(lsmean_cw)
#[1] 251   2

########-------CD
cd$Genotype <- as.character(cd$Genotype);
cd$p1 <- NA;
cd$p2 <- NA;
for (i in 1:nrow(cd)){
  tem <- sort(unlist(strsplit(cd$Genotype[i], split=" x ")));
  cd$p1[i] <- tem[1];
  cd$p2[i] <- tem[2];
}

cd <- cd[order(cd$p1, cd$p2),]
cd$ped <- paste(cd$p1, cd$p2, sep="|")
lsmean_cd <- ddply(cd, .(ped), summarise,
                   CD=mean(CD))
dim(lsmean_cd)
#[1] 251   2
######################################################################
nathan <- read.csv("Nathan_p1_subset.csv")
nam <- read.table("nam_parents", header=TRUE)
p <- c("B73", as.character(nam$parent))
p <- toupper(p);
nathan$INBRED <- toupper(nathan$INBRED)

pcob <- subset(nathan, INBRED %in% p)
pcob[pcob=="."] <- NA;
pcob$CW_Inbred <- as.numeric(as.character(pcob$CW_Inbred))
pcob$CL_Inbred <- as.numeric(as.character(pcob$CL_Inbred))
pcob$CD_Inbred <- as.numeric(as.character(pcob$CD_Inbred))
pcobmean <- ddply(pcob, .(INBRED), summarise,
              CW=mean(CW_Inbred, na.rm=T),
              CL=mean(CL_Inbred, na.rm=T),
              CD=mean(CD_Inbred, na.rm=T)
              )
b73 <- data.frame(INBRED="B73", CW=15.9, CL=115.8867, CD=25.59333)
pcobmean <- rbind(b73, pcobmean)
########################################################
#------------------ CL F1 and parents lsmean:
idx1 <- grep("@", lsmean_cl$ped)
f1cl <- lsmean_cl[-idx1,]

pcl <- lsmean_cl[idx1,]
pcl$ped <- gsub('@\\|NA',"", pcl$ped)

pcl$ped <- toupper(pcl$ped);
pcl2 <- merge(pcl, pcobmean[, c("INBRED", "CL")], by.x="ped", by.y="INBRED", all=T)

cor(subset(pcl2, !is.na(CL.x) & !is.na(CL.y))[,2:3])
#0.578

mymean <- function(x) mean(x, na.rm=T);
pcl2$CL <- apply(pcl2[,2:3], 1, mymean)
pcl2 <- pcl2[, c(1,4)]
names(pcl2)[1] <- "ped"

master_cl <- rbind(f1cl, pcl2)

#------------------ CW F1 and parents lsmean:
idx2 <- grep("@", lsmean_cw$ped)
f1cw <- lsmean_cw[-idx2,]

pcw <- lsmean_cw[idx2,]
pcw$ped <- gsub('@\\|NA',"", pcw$ped)

pcw$ped <- toupper(pcw$ped);
pcw2 <- merge(pcw, pcobmean[, c("INBRED", "CW")], by.x="ped", by.y="INBRED", all=T)

cor(subset(pcw2, !is.na(CW.x) & !is.na(CW.y))[,2:3])
#0.45

mymean <- function(x) mean(x, na.rm=T);
pcw2$CW <- apply(pcw2[,2:3], 1, mymean)
pcw2 <- pcw2[, c(1,4)]
names(pcw2)[1] <- "ped"

master_cw <- rbind(f1cw, pcw2)

#------------------ CD F1 and parents lsmean:
idx3 <- grep("@", lsmean_cd$ped)
f1cd <- lsmean_cd[-idx3,]

pcd <- lsmean_cd[idx3,]
pcd$ped <- gsub('@\\|NA',"", pcd$ped)

pcd$ped <- toupper(pcd$ped);
pcd2 <- merge(pcd, pcobmean[, c("INBRED", "CD")], by.x="ped", by.y="INBRED", all=T)

cor(subset(pcd2, !is.na(CD.x) & !is.na(CD.y))[,2:3])
#0.64

mymean <- function(x) mean(x, na.rm=T);
pcd2$CD <- apply(pcd2[,2:3], 1, mymean)
pcd2 <- pcd2[, c(1,4)]
names(pcd2)[1] <- "ped"

master_cd <- rbind(f1cd, pcd2)
###########################################
#### Merge the data for output
master_cob <- merge(master_cl, master_cw, by="ped", all=T)
master_cob <- merge(master_cob, master_cd, by="ped", all=T)

write.table(master_cob, "diallel_lsmean_cob_020112.csv", sep=",", row.names=FALSE, quote=FALSE)
## plot and output lsmean
idx2 <- grep("\\|", master_cd$ped);
f1 <- master_cd[idx2,]
p <- master_cd[-idx2, ]
plot(density(f1$CD), main="CD of Diallel and parents", xlab="KRN", lwd=3, col="blue")
lines(density(p$CD), lwd=3, lty=2, col="red")

