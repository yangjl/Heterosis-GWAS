### Jinliang Yang
### Updated: 6/13/2012
### using the IBM map to calculate the breakpoint

setwd("~/Documents/GWAS2_KRN/SNP/merged")
###### BxNAM RILs ###################
pheno <- read.table("~/Documents/GWAS2_KRN/pheno/krn_GenSel_6323_053112.txt", header=T)
bxril <- pheno[pheno$pop.=="BxRIL",]
bxril$Genotype <- gsub("B73x", "", bxril$Genotype)


########----------------Z017 data-----------------------------####################

z017 <- read.table("~/Documents/GWAS2_KRN/SNP/merged/tsnp_4impute_AGPv2.txt", header=TRUE, nrow=5)
z017id <- grep("Z017", names(z017), value=TRUE)

leaf <- read.csv("~/Documents/VariationDB/pheno/IBM_leaf_traits_20101025_KSU.csv", header=T)
leaf$MO.Number <- gsub("O", "0", leaf$MO.Number)

bxril <- merge(bxril, leaf[,1:2], by.x="Genotype", by.y="MO.Number", all.x=T)
tem1 <- bxril[is.na(bxril$Geno_Code),]
tem1$Geno_Code <- tem1$Genotype
tem2 <- bxril[!is.na(bxril$Geno_Code),]

bxril <- rbind(tem1, tem2)
bxril$Geno_Code <- gsub("M", "Z017EM", bxril$Geno_Code)
mxril <- subset(bxril, FID.=="Z017")

###########################################################
pullbx <- function(bpfile1="tagging_gbs50k_chr1.txt", bpfile2="IBM_impute/tagging_IBM136_chr1.txt", chr=1){
  tag1 <- read.table(bpfile1, header=TRUE);
  tag2 <- read.table(bpfile2, header=TRUE);
  
  
  tmp1 <- subset(tag1, id %in% bxril$Geno_Code)[,1:6];
  tmp2 <- subset(tag2, id %in% bxril$Geno_Code);
  
  tmp3 <- subset(tag1, id %in% mxril$Geno_Code)[,1:6];
  tmp4 <- subset(tag2, id %in% mxril$Geno_Code);
  
  out <- rbind(tmp1, tmp2)
  write.table(out, paste("BMxRILs/tagging_Bx692_chr", chr, ".txt", sep=""), row.names=FALSE, quote=FALSE, sep="\t")
  
  out2 <- rbind(tmp3, tmp4)
  write.table(out2, paste("BMxRILs/tagging_Mx289_chr", chr, ".txt", sep=""), row.names=FALSE, quote=FALSE, sep="\t")
  
  return( length(unique(out2$id)) )
}

pullbx(bpfile1="tagging_gbs50k_chr1.txt", bpfile2="IBM_impute/tagging_IBM136_chr1.txt", chr=1)
pullbx(bpfile1="tagging_gbs50k_chr2.txt", bpfile2="IBM_impute/tagging_IBM136_chr2.txt", chr=2)
pullbx(bpfile1="tagging_gbs50k_chr3.txt", bpfile2="IBM_impute/tagging_IBM136_chr3.txt", chr=3)
pullbx(bpfile1="tagging_gbs50k_chr4.txt", bpfile2="IBM_impute/tagging_IBM136_chr4.txt", chr=4)
pullbx(bpfile1="tagging_gbs50k_chr5.txt", bpfile2="IBM_impute/tagging_IBM136_chr5.txt", chr=5)
pullbx(bpfile1="tagging_gbs50k_chr6.txt", bpfile2="IBM_impute/tagging_IBM136_chr6.txt", chr=6)
pullbx(bpfile1="tagging_gbs50k_chr7.txt", bpfile2="IBM_impute/tagging_IBM136_chr7.txt", chr=7)
pullbx(bpfile1="tagging_gbs50k_chr8.txt", bpfile2="IBM_impute/tagging_IBM136_chr8.txt", chr=8)
pullbx(bpfile1="tagging_gbs50k_chr9.txt", bpfile2="IBM_impute/tagging_IBM136_chr9.txt", chr=9)
pullbx(bpfile1="tagging_gbs50k_chr10.txt", bpfile2="IBM_impute/tagging_IBM136_chr10.txt", chr=10)

##########   ----Bx -------####################################################################
./myimpute_beta4.pl --dsnp="clean_merged_chr1.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr1.txt" --cross="B73" --header=1 --mode=1 > BMxRILs/geno_Bx692_chr1.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr2.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr2.txt" --cross="B73" --header=1 --mode=1 > BMxRILs/geno_Bx692_chr2.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr3.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr3.txt" --cross="B73" --header=1 --mode=1 > BMxRILs/geno_Bx692_chr3.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr4.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr4.txt" --cross="B73" --header=1 --mode=1 > BMxRILs/geno_Bx692_chr4.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr5.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr5.txt" --cross="B73" --header=1 --mode=1 > BMxRILs/geno_Bx692_chr5.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr6.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr6.txt" --cross="B73" --header=1 --mode=1 > BMxRILs/geno_Bx692_chr6.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr7.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr7.txt" --cross="B73" --header=1 --mode=1 > BMxRILs/geno_Bx692_chr7.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr8.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr8.txt" --cross="B73" --header=1 --mode=1 > BMxRILs/geno_Bx692_chr8.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr9.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr9.txt" --cross="B73" --header=1 --mode=1 > BMxRILs/geno_Bx692_chr9.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr10.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr10.txt" --cross="B73" --header=1 --mode=1 > BMxRILs/geno_Bx692_chr10.GenSel

##########   ----Mx -------####################################################################
./myimpute_beta4.pl --dsnp="clean_merged_chr1.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr1.txt" --cross="Mo17" --header=1 --mode=1 > BMxRILs/geno_Mx289_chr1.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr2.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr2.txt" --cross="Mo17" --header=1 --mode=1 > BMxRILs/geno_Mx289_chr2.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr3.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr3.txt" --cross="Mo17" --header=1 --mode=1 > BMxRILs/geno_Mx289_chr3.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr4.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr4.txt" --cross="Mo17" --header=1 --mode=1 > BMxRILs/geno_Mx289_chr4.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr5.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr5.txt" --cross="Mo17" --header=1 --mode=1 > BMxRILs/geno_Mx289_chr5.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr6.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr6.txt" --cross="Mo17" --header=1 --mode=1 > BMxRILs/geno_Mx289_chr6.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr7.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr7.txt" --cross="Mo17" --header=1 --mode=1 > BMxRILs/geno_Mx289_chr7.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr8.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr8.txt" --cross="Mo17" --header=1 --mode=1 > BMxRILs/geno_Mx289_chr8.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr9.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr9.txt" --cross="Mo17" --header=1 --mode=1 > BMxRILs/geno_Mx289_chr9.GenSel
./myimpute_beta4.pl --dsnp="clean_merged_chr10.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr10.txt" --cross="Mo17" --header=1 --mode=1 > BMxRILs/geno_Mx289_chr10.GenSel






