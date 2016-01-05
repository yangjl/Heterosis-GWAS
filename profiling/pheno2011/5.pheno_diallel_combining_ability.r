# Jinliang Yang
# purpose: calculate the GCA and SCA
# updated: 3.7.2012

#Row        Barcode KRN Note.x         Pedigree   Genotype  Note.y
#3133  8757  11-8757-23 OP  16        10g-1028-21/2048 Tzi8 x P39 Diallel
# input data=d should contain trait and Genotype (Tzi x P39)

getca2 <- function(data=d, trait="KRN", square=FALSE){
  ### check the required cols
  if (!("Genotype" %in% names(data) & trait %in% names(data) )) 
    stop(paste("Make sure your data frame contains columns Genotype and ", trait, sep=""));
  data$Genotype <- as.character(data$Genotype)
  ### remove the self genotype
  idx <- grep("@", data$Genotype);
  if(sum(idx)>0) data <- data[-idx,];
  
  ### get the p1, p2 and reciprocal col
  if(!("p1" %in% names(data) & "p2" %in% names(data) & "reciprocal" %in% names(data) )){
    data$p1 <- NA;
    data$p2 <- NA;
    data$reciprocal <- NA;
    data$Genotype <- as.character(data$Genotype);
    for(i in 1:nrow(data)){
      tem <- unlist(strsplit(data$Genotype[i], split="x"));
      tem <- sort(tem)
      data$p1[i] <- tem[1];
      data$p2[i] <- tem[2];
      data$reciprocal[i] <- paste(tem[1], tem[2], sep="_")
    }
  }
  
  ### start the calculation according to Method 4, Model I
  uhat <- mean(data[,trait], na.rm=T);
  inbred <- sort(unique(c(data$p1, data$p2)));
  
  ### General Combining Ability
  gca <- data.frame(Genotype=inbred, GCA=NA);
  for(j in 1:length(inbred)){
      inbredj = inbred[j];
      ghat <- mean(subset(data, p1==inbredj | p2==inbredj)[, trait], na.rm=T)-uhat;
      gca[gca$Genotype==inbredj,]$GCA <- ghat;
  }
  
  ### Specific Combining Ability
  sca <- matrix(NA, nrow=length(inbred), ncol=length(inbred)+1);
  sca <- as.data.frame(sca)
  names(sca) <- c(trait, inbred)
  sca[,1] <- inbred
  
  specific <- unique(data$reciprocal);
  
  scatable <- data.frame()
  for(k in 1:length(specific)){
      tem2 <- unlist(strsplit(specific[k], split="_"));
      myp1 <- tem2[1];
      myp2 <- tem2[2];
      
      pheno <- subset(data, (p1==myp1 & p2==myp2) | (p1==myp2 & p2==myp1));
      shat <- mean(pheno[, trait], na.rm=T) - subset(gca, Genotype==myp1)$GCA - subset(gca, Genotype==myp2)$GCA - uhat;
    if(square){
      sca[sca[,1]==myp1, myp2] <- shat;
    } else {
      temid <- paste(myp1, myp2, sep="x");
      temsca <- data.frame(Genotype=temid, SCA=shat)
      scatable <- rbind(scatable, temsca)
    }
    
  }
  
  outputlist <- list();
  names(gca)[2] <- trait;
  outputlist[[1]] <- gca;
  names(scatable)[2] <- trait;
  if(square){
    outputlist[[2]] <- sca;
  } else{
    outputlist[[2]] <- scatable;
  }
  
  return(outputlist)
}
###############################################
mydata <- read.csv("pheno_diallel_master_BLUE.csv")

krnca <- getca2(data=mydata, trait="KRN", square=FALSE)
akwca <- getca2(data=mydata, trait="AKW", square=FALSE)
tkwca <- getca2(data=mydata, trait="TKW", square=FALSE)
kcca <- getca2(data=mydata, trait="KC", square=FALSE)

cdca <- getca2(data=mydata, trait="CD", square=FALSE)
cwca <- getca2(data=mydata, trait="CW", square=FALSE)
clca <- getca2(data=mydata, trait="CL", square=FALSE)

##########  GCA #################
gca_master <- merge(krnca[[1]], akwca[[1]], by="Genotype")
gca_master <- merge(gca_master, tkwca[[1]], by="Genotype")
gca_master <- merge(gca_master, kcca[[1]], by="Genotype")
gca_master <- merge(gca_master, cdca[[1]], by="Genotype")
gca_master <- merge(gca_master, cwca[[1]], by="Genotype")
gca_master <- merge(gca_master, clca[[1]], by="Genotype")

write.table(gca_master, "pheno_diallel_master_BLUE_GCA.csv", sep=",", row.names=FALSE, quote=FALSE)

##########  SCA #################
sca_master <- merge(krnca[[2]], akwca[[2]], by="Genotype")
sca_master <- merge(sca_master, tkwca[[2]], by="Genotype")
sca_master <- merge(sca_master, kcca[[2]], by="Genotype")
sca_master <- merge(sca_master, cdca[[2]], by="Genotype")
sca_master <- merge(sca_master, cwca[[2]], by="Genotype")
sca_master <- merge(sca_master, clca[[2]], by="Genotype")

write.table(sca_master, "pheno_diallel_master_BLUE_SCA.csv", sep=",", row.names=FALSE, quote=FALSE)





pairs(d[, c(4,7,9, 11:13)], text.panel = diag, upper.panel=panel.smooth, 
lower.panel=panel.cor, gap=0, main="KRN Correlation Plots of three obs")
  





#################################################################################
# make the design matrix
#####################################################################################################################################################
design <- read.table("/Users/yangjl/Documents/Heterosis_GWAS/SNP/diallel_4impute.txt", header=TRUE)

dm <- matrix(0, nrow=378, ncol=28+378)
dm <- as.data.frame(dm)
p <- as.character(unique(design$p1));
names(dm) <- c("Geno", p, as.character(design$geno))
dm$Geno <- as.character(design$geno);

for(i in 1:nrow(dm)){
  tem <- unlist(strsplit(dm$Geno[i], split="x"));
  p1 <- tem[1];
  p2 <- tem[2];
  dm[i, p1] <- 1;
  dm[i, p2] <- 1;
  dm[i, dm$Geno[i]] <- 1;
}

#################
trait <- read.csv("pheno_diallel_master_raw_040612.csv")
trait$Curtiss <- 0; # Only B73
trait$Dairy <- 0;
trait$Johnson <- 0;
trait$zumwalt <- 0;

for (i in 1:nrow(trait)){
  farm <- as.character(trait$Farm[i]);
  trait[i, farm] <- 1;
}

#############################
#FUNCTION to GCA and SCA
getca1 <- function(traitnm="TKW"){
  tem <- merge(trait[, c("Genotype2", traitnm, "Dairy", "Johnson", "zumwalt")], dm, by.x="Genotype2", by.y="Geno")
  tem <- tem[!is.na(tem[,2]),]
  yield <- tem[,2]
  
  idx1 <- numeric(ncol(tem))
  for(i in 3:ncol(tem)){
    idx1[i] <- sum(tem[,i])
  }
  idx <- which(idx1==0)
  
  dmatrix1 <- as.matrix(tem[, -idx]);
  fit1 <- lm(yield ~ dmatrix1)
  ped.hat1 <- fit1$coef
  ped.hat1 <- ped.hat1[!is.na(ped.hat1)]
  
  names(ped.hat1) <- gsub("dmatrix1", "", names(ped.hat1));
  names(ped.hat1)[1] <- "yhat"
  
  gcaout <- data.frame(Genotype=names(ped.hat1[1:29]), GCA=ped.hat1[1:29]);
  names(gcaout)[2] <- traitnm;
  
  scaout <- data.frame(Genotype=names(ped.hat1[30:227]), GCA=ped.hat1[30:227]);
  names(gcaout)[2] <- traitnm;
  
  
  output <- list(gca=gcaout, sca=scaout)
  return(output)
}

#############
krn1 <- getca1("KRN")

krntem <- merge(krnca[[2]], krn1[[2]], by="Genotype", all.x=TRUE)
krntem[is.na(krntem)] <- 0



kc <- getca("KC")
akw <- getca("AKW")
tkw <- getca("TKW")
cd <- getca("CD")
cw <- getca("CW")
cl <- getca("CL")

################################

getca(data=sample, trait="Yield")
#########################################################################################################################
x1 <- c("P1", "P1", "P1", "P1", "P2", "P2");
x2 <- c("P2", "P2", "P3", "P3", "P3", "P3")
yield <- c(10, 12, 20, 22, 14, 16)

sample <- data.frame(Parent1=x1, Parent2=x2, Yield=yield)

#####################################################################################################################################################
X <- matrix(c(1, 1, 0, 1, 0, 0,
                   1, 1, 0, 1, 0, 0,
                   1, 0, 1, 0, 1, 0,
                   1, 0, 1, 0, 1, 0,
                   0, 1, 1, 0, 0, 1,
                   0, 1, 1, 0, 0, 1), ncol=6, byrow=TRUE)

y <- sample$Yield
XX=t(X)%*%X 
library(MASS) 
XXgi=ginv(XX) 
Px=X%*%XXgi%*%t(X)
yhat=Px%*%y

fit <- lm(y~X)


y=c(4,1,3,5,3,3,1) 
X=matrix(c( 
  1,1,0,0,0,1,0,0, 
  1,1,0,0,0,0,1,0, 
  1,0,1,0,0,0,1,0, 
  1,0,1,0,0,0,0,1, 
  1,0,0,1,0,0,0,1, 
  1,0,0,0,1,1,0,0, 
  1,0,0,0,1,0,1,0 
  ),byrow=T,nrow=7) 
fit <- lm(y ~X)
XX=t(X)%*%X 
library(MASS) 
XXgi=ginv(XX) 
Px=X%*%XXgi%*%t(X)
yhat=Px%*%y

#################
data1 <- read.csv("~/Documents/Heterosis_GWAS/pheno2011/LL_KN_AKW_test.csv");

###########################################################
### Delete duplicated input and missing data
file1 <- read.csv("LL_KRN_test.csv")
dim(file1)
#8678 8

diallel <- subset(file1, Note.y=="Diallel" | Note.y=="NAM Filler");
dim(diallel)
#[1] 2850    8
#summary(diallel)

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


data2 <- read.csv("~/Documents/Heterosis_GWAS/pheno2011/LL_KRN_test.csv");
data3 <- read.csv("~/Documents/Heterosis_GWAS/pheno2011/LL_cob.test.csv");

data1 <- subset(data1, )

dm2 <- merge(sample, dm, by.x="Genotype", by.y="Geno")

yield <- dm2$Yield
dmatrix <- as.matrix(dm2[, 7:52]);

fit <- lm(yield ~ dmatrix)




