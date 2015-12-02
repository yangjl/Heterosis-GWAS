# Jinliang Yang
# update: 6.5.2012

### SNP extension ###
setwd("/Users/yangjl/Documents/GWAS2_KRN/SNP/merged")
source("bpExtension.R") 
bpExtension(bpfile="bp_tagging1000.txt", GBSfile="tsnp_GBS50k_AGPv2.txt", parellel=TRUE, range=1:600, temfile="temfile1", outputfile="bp_GBS50k_1.txt")
### 4892 rils
### SNP extension ###
setwd("/Users/yangjl/Documents/GWAS2_KRN/SNP/merged")
source("bpExtension.R") 
bpExtension(bpfile="bp_tagging1000.txt", GBSfile="tsnp_GBS50k_AGPv2.txt", parellel=TRUE, range=601:1200, temfile="temfile2", outputfile="bp_GBS50k_2.txt")
### 4892 rils
setwd("/Users/yangjl/Documents/GWAS2_KRN/SNP/merged")
source("bpExtension.R") 
bpExtension(bpfile="bp_tagging1000.txt", GBSfile="tsnp_GBS50k_AGPv2.txt", parellel=TRUE, range=1201:1800, temfile="temfile3", outputfile="bp_GBS50k_3.txt")
### 4892 rils
setwd("/Users/yangjl/Documents/GWAS2_KRN/SNP/merged")
source("bpExtension.R") 
bpExtension(bpfile="bp_tagging1000.txt", GBSfile="tsnp_GBS50k_AGPv2.txt", parellel=TRUE, range=1801:2400, temfile="temfile4", outputfile="bp_GBS50k_4.txt")
### 4892 rils
setwd("/Users/yangjl/Documents/GWAS2_KRN/SNP/merged")
source("bpExtension.R") 
bpExtension(bpfile="bp_tagging1000.txt", GBSfile="tsnp_GBS50k_AGPv2.txt", parellel=TRUE, range=2401:3000, temfile="temfile5", outputfile="bp_GBS50k_5.txt")
### 4892 rils
setwd("/Users/yangjl/Documents/GWAS2_KRN/SNP/merged")
source("bpExtension.R") 
bpExtension(bpfile="bp_tagging1000.txt", GBSfile="tsnp_GBS50k_AGPv2.txt", parellel=TRUE, range=3001:3600, temfile="temfile6", outputfile="bp_GBS50k_6.txt")
### 4892 rils
setwd("/Users/yangjl/Documents/GWAS2_KRN/SNP/merged")
source("bpExtension.R") 
bpExtension(bpfile="bp_tagging1000.txt", GBSfile="tsnp_GBS50k_AGPv2.txt", parellel=TRUE, range=3601:4200, temfile="temfile7", outputfile="bp_GBS50k_7.txt")
### 4892 rils
setwd("/Users/yangjl/Documents/GWAS2_KRN/SNP/merged")
source("bpExtension.R") 
bpExtension(bpfile="bp_tagging1000.txt", GBSfile="tsnp_GBS50k_AGPv2.txt", parellel=TRUE, range=4201:4892, temfile="temfile8", outputfile="bp_GBS50k_8.txt")
### 4892 rils

###########
gbs50k1 <- read.csv("bp_GBS50k_1.txt", header=F)
gbs50k2 <- read.csv("bp_GBS50k_2.txt", header=F)
gbs50k3 <- read.csv("bp_GBS50k_3.txt", header=F)
gbs50k4 <- read.csv("bp_GBS50k_4.txt", header=F)
gbs50k5 <- read.csv("bp_GBS50k_5.txt", header=F)
gbs50k6 <- read.csv("bp_GBS50k_6.txt", header=F)
gbs50k7 <- read.csv("bp_GBS50k_7.txt", header=F)
gbs50k8 <- read.csv("bp_GBS50k_8.txt", header=F)

gbs50k <- rbind(gbs50k1, gbs50k2,gbs50k3, gbs50k4,gbs50k5, gbs50k6,gbs50k7, gbs50k8);
names(gbs50k) <- c("id", "chr", "pos1", "pos2", "snp1", "snp2", "dis")
#gbs02 <- subset(gbs50k, snp1 != snp2)
length(unique(gbs50k$id)) #4892

chr1 <- subset(gbs50k, chr==1)
chr2 <- subset(gbs50k, chr==2)
chr3 <- subset(gbs50k, chr==3)
chr4 <- subset(gbs50k, chr==4)
chr5 <- subset(gbs50k, chr==5)
chr6 <- subset(gbs50k, chr==6)
chr7 <- subset(gbs50k, chr==7)
chr8 <- subset(gbs50k, chr==8)
chr9 <- subset(gbs50k, chr==9)
chr10 <- subset(gbs50k, chr==10)

write.table(chr1, "tagging_gbs50k_chr1.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr2, "tagging_gbs50k_chr2.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr3, "tagging_gbs50k_chr3.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr4, "tagging_gbs50k_chr4.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr5, "tagging_gbs50k_chr5.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr6, "tagging_gbs50k_chr6.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr7, "tagging_gbs50k_chr7.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr8, "tagging_gbs50k_chr8.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr9, "tagging_gbs50k_chr9.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr10, "tagging_gbs50k_chr10.txt", sep="\t", row.names=FALSE, quote=FALSE)













