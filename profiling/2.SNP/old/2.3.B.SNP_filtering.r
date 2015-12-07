# Merging the two set of SNP
# update: 04/07/2011
# updated: 10/10/2011

##################################################################
# PLINK Filtering by missing rate and MAF
##################################################################
setwd("/Users/yangjl/Documents/KRN_GWAS/SNP/merged")

snpfile <- read.delim("SNP_founders_012612.txt", header=TRUE)
mysnp <- snpfile
dim(mysnp)
#[1] 3380272      31

#################################################################
# get the fine scale map
#################################################################
#SNP Chr Physical
map <- mysnp[,1:3]
names(map) <- c("SNP", "Chr", "Physical")

write.table(map, "../finemap/merged_012612.txt", sep="\t", row.names=FALSE, quote=FALSE)

##################################################################
##### get rid of the additional alleles
Apex <- mysnp[mysnp$Set==12,]

Apex$N <- apply(Apex[, 4:30], 1, function(x) sum(x=="N")) 
Apex$A <- apply(Apex[, 4:30], 1, function(x) sum(x=="A")) 
Apex$T <- apply(Apex[, 4:30], 1, function(x) sum(x=="T")) 
Apex$C <- apply(Apex[, 4:30], 1, function(x) sum(x=="C")) 
Apex$G <- apply(Apex[, 4:30], 1, function(x) sum(x=="G"))
Apex$M <- apply(Apex[, 4:30], 1, function(x) sum(x=="-"))
Apex$P <- apply(Apex[, 4:30], 1, function(x) sum(x=="+"))

Apex$count0 <- apply(Apex[, 33:38], 1, function(x) sum(x==0))

Apex_4alle <- Apex[Apex$count0 ==1,]
Apex1 <- Apex[Apex$count0!=4,]
Apex2 <- Apex[Apex$count0==4,]

exclude <- as.character(Apex1$rs)
tem <- mysnp[!(mysnp$rs %in% exclude),]
mysnp <- tem[order(tem$chr, tem$chr),]
dim(mysnp)
#[1] 3314971      31


##################################################################


#-------->output the .tmap
nm <- names(mysnp)[4:30]
mysnp.tmap<- data.frame(family= nm, individual=c(1:27), paternal=0, maternal=0, sex=0, pheno=0)
write.table(mysnp.tmap, "merged_012612.tfam", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

#-------->output the .tped
mysnp <- mysnp[order(mysnp$chr, mysnp$pos),]
mysnp$Set <- 1:nrow(mysnp)
### Haploid to Diploid:
for (i in 4:30){
	mysnp[,i] <- paste(mysnp[,i], mysnp[,i], sep=" ")
}

mysnp$genpos <- mysnp$pos/10000
mysnp <- mysnp[, c(1:2,32,3:30)]
#chr
#rs
#genetic pos
#physical pos
#genotype...

write.table(mysnp, "merged_012612.tped", sep="\t", quote=FALSE, col.names=FALSE,row.names=FALSE)


################################################################################
# Then PLINK Command:

# 1. make the bed
plink --tped merged_012612.tped --tfam merged_012612.tfam --missing-genotype N --make-bed --out merged_012612


# 2. merge the file to one
plink --bfile merged_012612 --missing --out merged_012612
plink --bfile merged_012612 --freq --out merged_012612

################################################################################
# start the filtering based on allele missing rate and MAF
missing <- read.table("merged_012612.lmiss", header=TRUE)
dim(missing[missing$F_MISS <= 0.6,])
m1 <- missing[missing$F_MISS <= 0.6,]

# Because B73 allele can be enlarged in the following step, so, we can not filtering based on the parental MAF
freq <- read.table("merged_012612.frq", header=TRUE)
dim(freq[freq$MAF > 0.1,])
freq1 <- freq[freq$MAF > 0.1,]

################################################################################
#for imputation use
setwd("/Users/yangjl/Documents/KRN_GWAS/SNP/merged")

mysnp <- read.delim("SNP_founders_012612.txt", header=TRUE)
#map <- read.delim("/Users/yangjl/Documents/SNP/finemap/HapApex_m_040711.txt.gen",header=TRUE)
#mysnp <- merge(mysnp, map, by.x="rs", by.y="marker")

#mysnp <- mysnp[,c(1:2,34,3:30)]
#names(mysnp)[2] <- "chr"
#names(mysnp)[3] <- "genetpos"

fdsnp <- mysnp[mysnp$rs %in% freq1$SNP,] # used the missing rate <= 0.6 threshold and MAF > 0.1
fdsnp <- fdsnp[order(fdsnp$chr, fdsnp$pos),]
dim(fdsnp)
#[1]2059242      31

########### remove the duplicated SNPs
fdsnp$uid <- paste(fdsnp$chr, fdsnp$pos, sep="_");
fdsnp <- fdsnp[!duplicated(fdsnp$uid),]
dim(fdsnp)
#[1] 2058739      32

############################################
fdsnp.map <- fdsnp[,1:3]
fdsnp.map$genetpos <- fdsnp.map$pos/10000
fdsnp.map <- fdsnp.map[,c(1,2,4,3)]

#for imputation use
write.table(fdsnp[,-31:-32], "fdsnp_2M_012612.dsnp", sep="\t", quote=FALSE, row.names=FALSE)

#output the .map
write.table(fdsnp.map, "fdsnp_2M_012612_v1.map", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

tem2 <- grep("PS", fdsnp.map$rs)




