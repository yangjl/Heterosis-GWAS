### Jinliang Yang
### Updated: 6/13/2012
### using the IBM map to calculate the breakpoint


###### MAGI MAP of 1017 SNP for 292 IBM RILs ###################
# MAP
geno <- read.csv("/Users/yangjl/Documents/QTL_seg/QTL_1017SNP/geno.ibm.snp_1016.csv")
dim(geno)

geno <- t(geno)
colnm <- c("Chr", "GenetPos",geno[1,c(-1,-2)])
geno <- geno[-1,]
rownm <- row.names(geno)
geno <- apply(geno, 2, as.numeric)
geno <- as.data.frame(geno)
names(geno) <- colnm
geno <- cbind(data.frame(Marker=rownm), geno)

map <- read.table("~/Documents/QTL_seg/QTL_1017SNP/MAGI_SNP1016.imputed.txt", header=T)
magi <- merge(map[,c(1,3,6)], geno, by.x="marker", by.y="Marker")
magi <- magi[, c(1,2, 4, 5, 3, 6:295)]
magi <- magi[order(magi$Chr, magi$imputepos),]

### change code B73 0->2; Mo17 1-> 0; missing NA-> -9;
tem <- magi[, 6:295] 
tem[tem==0] <- 2;
tem[tem==1] <- 0;
tem[is.na(tem)] <- -9;

magi[, 6:295] <- tem;
magiid <- names(magi)

########----------------Z017 data-----------------------------####################

z017 <- read.table("~/Documents/GWAS2_KRN/SNP/merged/tsnp_4impute_AGPv2.txt", header=TRUE, nrow=5)
z017id <- grep("Z017", names(z017), value=TRUE)

leaf <- read.csv("~/Documents/VariationDB/pheno/IBM_leaf_traits_20101025_KSU.csv", header=T)
leaf$MO.Number <- gsub("O", "0", leaf$MO.Number)


#################################################################
list <- unique(c(leaf$MO.Number,  magiid[-1:-5]));
#329

###--- Genotyped by 50k GBS
idx1 <- which(leaf$MO.Number %in% list)
list1<- list[-idx1]

###########################################################
magi2 <- magi[, c(magiid[1:5], list1)]
names(magi2)[1:5] <- names(z017)[1:5]
magi2 <- magi2[order(magi2$ch, magi2$pos),]

write.table(magi2, "tagging_IBM_136_4impute.txt", sep="\t", row.names=FALSE, quote=FALSE)

######
./myimpute_calbp.pl --tsnp="tagging_IBM_136_4impute.txt" --snpwindow=0 --start=5 --end=0 > IBM_impute/tagging_IBM_MAGI.txt


bp_magi <- read.table("~/Documents/GWAS2_KRN/SNP/merged/IBM_impute/tagging_IBM_MAGI.txt", header=FALSE)

names(bp_magi) <- c("id", "chr", "pos1", "pos2", "snp1", "snp2")
#gbs02 <- subset(gbs50k, snp1 != snp2)
length(unique(bp_magi$id)) #136
bp_magi$id <- gsub("M", "Z017EM", bp_magi$id)

chr1 <- subset(bp_magi, chr==1)
chr2 <- subset(bp_magi, chr==2)
chr3 <- subset(bp_magi, chr==3)
chr4 <- subset(bp_magi, chr==4)
chr5 <- subset(bp_magi, chr==5)
chr6 <- subset(bp_magi, chr==6)
chr7 <- subset(bp_magi, chr==7)
chr8 <- subset(bp_magi, chr==8)
chr9 <- subset(bp_magi, chr==9)
chr10 <- subset(bp_magi, chr==10)

write.table(chr1, "tagging_IBM136_chr1.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr2, "tagging_IBM136_chr2.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr3, "tagging_IBM136_chr3.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr4, "tagging_IBM136_chr4.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr5, "tagging_IBM136_chr5.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr6, "tagging_IBM136_chr6.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr7, "tagging_IBM136_chr7.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr8, "tagging_IBM136_chr8.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr9, "tagging_IBM136_chr9.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(chr10, "tagging_IBM136_chr10.txt", sep="\t", row.names=FALSE, quote=FALSE)

######-------------RIL---------------------###########################################
./myimpute_beta4.pl --dsnp="clean_merged_chr1.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr1.txt" --cross="RIL" --header=1 --mode=1 > IBM_impute/geno_IBM136_chr1.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr2.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr2.txt" --cross="RIL" --header=1 --mode=1 > IBM_impute/geno_IBM136_chr2.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr3.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr3.txt" --cross="RIL" --header=1 --mode=1 > IBM_impute/geno_IBM136_chr3.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr4.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr4.txt" --cross="RIL" --header=1 --mode=1 > IBM_impute/geno_IBM136_chr4.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr5.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr5.txt" --cross="RIL" --header=1 --mode=1 > IBM_impute/geno_IBM136_chr5.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr6.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr6.txt" --cross="RIL" --header=1 --mode=1 > IBM_impute/geno_IBM136_chr6.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr7.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr7.txt" --cross="RIL" --header=1 --mode=1 > IBM_impute/geno_IBM136_chr7.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr8.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr8.txt" --cross="RIL" --header=1 --mode=1 > IBM_impute/geno_IBM136_chr8.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr9.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr9.txt" --cross="RIL" --header=1 --mode=1 > IBM_impute/geno_IBM136_chr9.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr10.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr10.txt" --cross="RIL" --header=1 --mode=1 > IBM_impute/geno_IBM136_chr10.GenSel

######-------------BxRIL---------------------###########################################






