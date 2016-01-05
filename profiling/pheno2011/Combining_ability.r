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
  outputlist[[1]] <- gca;
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
trait$u <- 1
getca1 <- function(traitnm="KRN"){
  tem <- merge(trait[, c("Genotype2", traitnm,"u", "Dairy", "Johnson", "zumwalt")], dm, by.x="Genotype2", by.y="Geno")
  tem <- tem[!is.na(tem[,2]),]
  yield <- tem[,2]
  dmatrix1 <- as.matrix(tem[, 3:5]);
  
  dmatrix1 <- apply(dmatrix1,2, as.factor)
  
  
  dmatrix1 <- as.matrix(tem[, 6:410]);
  fit2 <- lm(yield ~ dmatrix1)
  ped.hat1 <- fit1$coef
  ped.hat1 <- ped.hat1[!is.na(ped.hat1)]
  
  names(ped.hat1) <- gsub("dmatrix1", "", names(ped.hat1));
  names(ped.hat1)[1] <- "yhat"
  
  output <- list(gca=ped.hat1[1:29], sca=ped.hat1[30:227])
  return(output)
}

#############
krn1 <- getca1("KRN")
kc <- getca("KC")
akw <- getca("AKW")
tkw <- getca("TKW")
cd <- getca("CD")
cw <- getca("CW")
cl <- getca("CL")

################################

getca(data=sample, trait="Yield")
#########################################################################################################################
x1 <- c(rep("P1",8), rep("P2",7), rep("P3",6), rep("P4",5), rep("P5",4), rep("P6",3), rep("P7",2), rep("P8",1));
x2 <- c(paste("P",2:9,sep=""), paste("P",3:9,sep=""), paste("P",4:9,sep=""), paste("P",5:9,sep=""), paste("P",6:9,sep=""),
        paste("P",7:9,sep=""), paste("P",8:9,sep=""), paste("P",9,sep=""))
yield <- c(240, 260, 230.4, 257, 241.5, 266.9, 240.1, 300.4, 
               209, 217.3, 233.1, 229.5, 266.9, 216.3, 214.2,
                    183.7, 253.7, 250.1, 268.8, 222.3, 252.1,
                           233.8, 213.7, 255.7, 197.4, 281,
                                  206.8, 272.2, 242.9, 260.8,
                                         261.8, 270.3, 283.9,
                                                273.2, 302.2,
                                                       259.8)
cw <- c(31.8, 34.7, 32.3, 45, 39, 35.1, 35.7, 40.1,
              27.9, 30.8, 39.6, 33.1, 30.9, 30.9, 28.5,
                    25.2, 41.4, 35.5, 34.9, 32.1, 32.4,
                          42.6, 35.7, 35.6, 32.7, 41.3,
                                40.1, 43.6, 41.8, 44.2,
                                      39.1, 43.5, 41.5,
                                            38.3, 41.1,
                                                  35.2)
scw <- c(208.2, 225.3, 198.1, 212, 202.5, 231.8, 204.4, 260.3,
                181.1, 186.5, 193.5, 196.4, 236, 185.4, 185.7,
                       158.5, 212.3, 214.6, 233.9, 190.2, 219.7,
                              191.2, 178, 220.1, 164.7, 239.7,
                                     166.7, 228.6, 201.1, 216.6,
                                            222.7, 226.8, 242.4,
                                                   234.9, 261.1,
                                                          224.6)

sample <- data.frame(P1=x1, P2=x2, Yield=yield, CW=cw, SCW=scw)
sample$Genotype <- paste(sample$P1, sample$P2, sep=" x ")


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




