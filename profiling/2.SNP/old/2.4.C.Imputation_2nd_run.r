# Jinliang Yang
# purpose: Merged SNP Imputation
# Last update: 6.12.2012

########################################################################################################
# prepare the selected SNP for 2nd round imputation
##################################################################################################
# with PAV SNP
setwd("")
mf <- read.csv("~/Documents/GWAS2_KRN/Method/GenSel/KRN_072512/bayes_GenSel_072712.csv", header=TRUE);

tab5rows <- read.table("~/Documents/GWAS2_KRN/SNP/merged/cleaned/clean_merged_chr1.dsnp", header=TRUE, nrows=5)
classes <- sapply(tab5rows, class)

chr1 <- read.table("~/Documents/GWAS2_KRN/SNP/merged/cleaned/clean_merged_chr1.dsnp", header=TRUE, colClasses=classes)
chr2 <- read.table("~/Documents/GWAS2_KRN/SNP/merged/cleaned/clean_merged_chr2.dsnp", header=TRUE, colClasses=classes)
chr3 <- read.table("~/Documents/GWAS2_KRN/SNP/merged/cleaned/clean_merged_chr3.dsnp", header=TRUE, colClasses=classes)
chr4 <- read.table("~/Documents/GWAS2_KRN/SNP/merged/cleaned/clean_merged_chr4.dsnp", header=TRUE, colClasses=classes)
chr5 <- read.table("~/Documents/GWAS2_KRN/SNP/merged/cleaned/clean_merged_chr5.dsnp", header=TRUE, colClasses=classes)
chr6 <- read.table("~/Documents/GWAS2_KRN/SNP/merged/cleaned/clean_merged_chr6.dsnp", header=TRUE, colClasses=classes)
chr7 <- read.table("~/Documents/GWAS2_KRN/SNP/merged/cleaned/clean_merged_chr7.dsnp", header=TRUE, colClasses=classes)
chr8 <- read.table("~/Documents/GWAS2_KRN/SNP/merged/cleaned/clean_merged_chr8.dsnp", header=TRUE, colClasses=classes)
chr9 <- read.table("~/Documents/GWAS2_KRN/SNP/merged/cleaned/clean_merged_chr9.dsnp", header=TRUE, colClasses=classes)
chr10 <- read.table("~/Documents/GWAS2_KRN/SNP/merged/cleaned/clean_merged_chr10.dsnp", header=TRUE, colClasses=classes);

chrpav <- read.table("~/Documents/GWAS2_KRN/SNP/merged/all_merged/pav_SNP_all.dsnp", header=TRUE, colClasses=classes);
chrpav$rs <- paste("pav", chrpav$rs)
chrpav$rs <- gsub(" ", "", chrpav$rs)

chr <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chrpav);
dim(chr)
#[1] 14401282       31

mychr <- chr[chr$rs %in% mf$ID,]
dim(mychr)
#[1] 1440128      31

write.table(mychr, "~/Documents/GWAS2_KRN/SNP/merged/krn_mf_1.4M.dsnp", row.names=FALSE, quote=FALSE, sep="\t")
write.table(mychr[,1:4], "~/Documents/GWAS2_KRN/SNP/merged/geno_krn_mf_PLINK.map", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
##################################################################################################
# without PAV SNP
setwd("")
mf <- read.csv("~/Documents/GWAS2_KRN/Method/GenSel/KRN_072512/bayes_GenSel_072712_nopav.csv", header=TRUE);
mychr <- read.table("~/Documents/GWAS2_KRN/SNP/merged/krn_mf_1.4M.dsnp", header=T)
mf2 <- mychr[mychr$rs %in% mf$ID,];
dim(mf2)
#[1] 1296628      31
write.table(mf2, "~/Documents/GWAS2_KRN/SNP/merged/krn_mf_1.3M.dsnp", row.names=FALSE, quote=FALSE, sep="\t")
write.table(mf2[,1:4], "~/Documents/GWAS2_KRN/SNP/merged/geno_krn_mf_1.3M_PLINK.map", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


########################################################################################################
# IMPUTE for GenSel
##################################################################################################
./myimpute_beta4.pl --dsnp="krn_mf_1.4M.dsnp" --bpfile="cleaned/tagging_gbs50k_all.txt" --header=1 --mode=1 --cross="RIL" > geno_krn_mf_GenSel.txt
./myimpute_beta4.pl --dsnp="krn_mf_1.4M.dsnp" --bpfile="IBM_impute/tagging_IBM136_all.txt" --header=2 --mode=1 --cross="RIL" >> geno_krn_mf_GenSel.txt
./myimpute_beta4.pl --dsnp="krn_mf_1.4M.dsnp" --bpfile="BMxRILs/tagging_Bx692_all.txt" --header=2 --mode=1 --cross="B73" >> geno_krn_mf_GenSel.txt
./myimpute_beta4.pl --dsnp="krn_mf_1.4M.dsnp" --bpfile="BMxRILs/tagging_Mx289_all.txt" --header=2 --mode=1 --cross="Mo17" >> geno_krn_mf_GenSel.txt
./myimpute_diallel.pl --diallel="Diallel/diallel_4impute.txt" --dsnp="krn_mf_1.4M.dsnp" --header=2 --mode=1 >> geno_krn_mf_GenSel.txt

./newgen2bin -i geno_krn_mf_GenSel.txt

######################################################################################################
# IMPUTE for PLINK
##################################################################################################
./myimpute_beta4.pl --dsnp="krn_mf_1.4M.dsnp" --bpfile="cleaned/tagging_gbs50k_all.txt" --header=2 --mode=2 --cross="RIL" > geno_krn_mf_PLINK.ped
./myimpute_beta4.pl --dsnp="krn_mf_1.4M.dsnp" --bpfile="IBM_impute/tagging_IBM136_all.txt" --header=2 --mode=2 --cross="RIL" >> geno_krn_mf_PLINK.ped
./myimpute_beta4.pl --dsnp="krn_mf_1.4M.dsnp" --bpfile="BMxRILs/tagging_Bx692_all.txt" --header=2 --mode=2 --cross="B73" >> geno_krn_mf_PLINK.ped
./myimpute_beta4.pl --dsnp="krn_mf_1.4M.dsnp" --bpfile="BMxRILs/tagging_Mx289_all.txt" --header=2 --mode=2 --cross="Mo17" >> geno_krn_mf_PLINK.ped
./myimpute_diallel.pl --dsnp="krn_mf_1.4M.dsnp" --diallel="Diallel/diallel_4impute.txt" --header=2 --mode=2 >> geno_krn_mf_PLINK.ped









