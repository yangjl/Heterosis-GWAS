# Jinliang Yang
# purpose: Merged SNP Imputation
# Last update: 6.12.2012

########################################################################################################
# NAM RILs Imputation
##################################################################################################
setwd("~/Documents/GWAS2_KRN/SNP/merged/cleaned")
tped5rows <- read.table("cleaned_SNP_chr1.tped", header=F, nrow=5)
classes <- sapply(tped5rows, class)
nm <- c("chr", "rs","allele", "pos","B73", "Z001","Z002", "Z003","Z004", 
  "Z005","Z006","Z007","Z008","Z009","Z010",
  "Z011","Z012","Z013","Z014","Z015","Z016",
  "Z017","Z018","Z019","Z020","Z021","Z022",
  "Z023","Z024","Z025","Z026")

tped2dsnp <- function(tpedfile="cleaned_SNP_chr1.tped", outfile="clean_merged_chr1.dsnp"){
  tped1 <- read.table(tpedfile, header=FALSE, colClasses = classes)
  tped1 <- tped1[, c(1:4, 5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57)]
  names(tped1) <- nm
  tem <- tped1[, 5:31]
  tem <- apply(tem, 2, as.character)
  tem[tem==0] <- "N"
  tped1[,5:31] <- tem
  write.table(tped1, outfile, sep="\t", row.names=FALSE, quote=FALSE)
  return(dim(tped1))
}
  
tped2dsnp(tpedfile="cleaned_SNP_chr1.tped", outfile="clean_merged_chr1.dsnp")
tped2dsnp(tpedfile="cleaned_SNP_chr2.tped", outfile="clean_merged_chr2.dsnp")
tped2dsnp(tpedfile="cleaned_SNP_chr3.tped", outfile="clean_merged_chr3.dsnp")
tped2dsnp(tpedfile="cleaned_SNP_chr4.tped", outfile="clean_merged_chr4.dsnp")
tped2dsnp(tpedfile="cleaned_SNP_chr5.tped", outfile="clean_merged_chr5.dsnp")
tped2dsnp(tpedfile="cleaned_SNP_chr6.tped", outfile="clean_merged_chr6.dsnp")
tped2dsnp(tpedfile="cleaned_SNP_chr7.tped", outfile="clean_merged_chr7.dsnp")
tped2dsnp(tpedfile="cleaned_SNP_chr8.tped", outfile="clean_merged_chr8.dsnp")
tped2dsnp(tpedfile="cleaned_SNP_chr9.tped", outfile="clean_merged_chr9.dsnp")
tped2dsnp(tpedfile="cleaned_SNP_chr10.tped", outfile="clean_merged_chr10.dsnp")

###################################################################################
# NAM RILs
###################################################################################
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr1.dsnp" --bpfile="cleaned/tagging_gbs50k_chr1.txt" --header=2 --mode=2 --cross="RIL" > geno_chr1.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr2.dsnp" --bpfile="cleaned/tagging_gbs50k_chr2.txt" --header=2 --mode=2 --cross="RIL" > geno_chr2.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr3.dsnp" --bpfile="cleaned/tagging_gbs50k_chr3.txt" --header=2 --mode=2 --cross="RIL" > geno_chr3.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr4.dsnp" --bpfile="cleaned/tagging_gbs50k_chr4.txt" --header=2 --mode=2 --cross="RIL" > geno_chr4.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr5.dsnp" --bpfile="cleaned/tagging_gbs50k_chr5.txt" --header=2 --mode=2 --cross="RIL" > geno_chr5.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr6.dsnp" --bpfile="cleaned/tagging_gbs50k_chr6.txt" --header=2 --mode=2 --cross="RIL" > geno_chr6.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr7.dsnp" --bpfile="cleaned/tagging_gbs50k_chr7.txt" --header=2 --mode=2 --cross="RIL" > geno_chr7.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr8.dsnp" --bpfile="cleaned/tagging_gbs50k_chr8.txt" --header=2 --mode=2 --cross="RIL" > geno_chr8.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr9.dsnp" --bpfile="cleaned/tagging_gbs50k_chr9.txt" --header=2 --mode=2 --cross="RIL" > geno_chr9.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr10.dsnp" --bpfile="cleaned/tagging_gbs50k_chr10.txt" --header=2 --mode=2 --cross="RIL" > geno_chr10.PLINK

###################################################################################
# IBM RILs
###################################################################################
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr1.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr1.txt" --header=2 --mode=2 --cross="RIL" > IBM_impute/geno_PLINK_ibm136_chr1.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr2.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr2.txt" --header=2 --mode=2 --cross="RIL" > IBM_impute/geno_PLINK_ibm136_chr2.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr3.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr3.txt" --header=2 --mode=2 --cross="RIL" > IBM_impute/geno_PLINK_ibm136_chr3.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr4.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr4.txt" --header=2 --mode=2 --cross="RIL" > IBM_impute/geno_PLINK_ibm136_chr4.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr5.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr5.txt" --header=2 --mode=2 --cross="RIL" > IBM_impute/geno_PLINK_ibm136_chr5.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr6.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr6.txt" --header=2 --mode=2 --cross="RIL" > IBM_impute/geno_PLINK_ibm136_chr6.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr7.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr7.txt" --header=2 --mode=2 --cross="RIL" > IBM_impute/geno_PLINK_ibm136_chr7.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr8.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr8.txt" --header=2 --mode=2 --cross="RIL" > IBM_impute/geno_PLINK_ibm136_chr8.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr9.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr9.txt" --header=2 --mode=2 --cross="RIL" > IBM_impute/geno_PLINK_ibm136_chr9.PLINK
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr10.dsnp" --bpfile="IBM_impute/tagging_IBM136_chr10.txt" --header=2 --mode=2 --cross="RIL" > IBM_impute/geno_PLINK_ibm136_chr10.PLINK

###################################################################################
# B73 x NAM RILs
###################################################################################
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr1.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr1.txt" --header=2 --mode=2 --cross="B73" > BMxRILs/geno_PLINK_Bx_chr1.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr2.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr2.txt" --header=2 --mode=2 --cross="B73" > BMxRILs/geno_PLINK_Bx_chr2.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr3.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr3.txt" --header=2 --mode=2 --cross="B73" > BMxRILs/geno_PLINK_Bx_chr3.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr4.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr4.txt" --header=2 --mode=2 --cross="B73" > BMxRILs/geno_PLINK_Bx_chr4.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr5.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr5.txt" --header=2 --mode=2 --cross="B73" > BMxRILs/geno_PLINK_Bx_chr5.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr6.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr6.txt" --header=2 --mode=2 --cross="B73" > BMxRILs/geno_PLINK_Bx_chr6.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr7.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr7.txt" --header=2 --mode=2 --cross="B73" > BMxRILs/geno_PLINK_Bx_chr7.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr8.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr8.txt" --header=2 --mode=2 --cross="B73" > BMxRILs/geno_PLINK_Bx_chr8.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr9.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr9.txt" --header=2 --mode=2 --cross="B73" > BMxRILs/geno_PLINK_Bx_chr9.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr10.dsnp" --bpfile="BMxRILs/tagging_Bx692_chr10.txt" --header=2 --mode=2 --cross="B73" > BMxRILs/geno_PLINK_Bx_chr10.ped 

###################################################################################
# Mo17 x NAM RILs
###################################################################################
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr1.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr1.txt" --header=2 --mode=2 --cross="Mo17" > BMxRILs/geno_PLINK_Mx_chr1.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr2.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr2.txt" --header=2 --mode=2 --cross="Mo17" > BMxRILs/geno_PLINK_Mx_chr2.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr3.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr3.txt" --header=2 --mode=2 --cross="Mo17" > BMxRILs/geno_PLINK_Mx_chr3.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr4.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr4.txt" --header=2 --mode=2 --cross="Mo17" > BMxRILs/geno_PLINK_Mx_chr4.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr5.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr5.txt" --header=2 --mode=2 --cross="Mo17" > BMxRILs/geno_PLINK_Mx_chr5.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr6.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr6.txt" --header=2 --mode=2 --cross="Mo17" > BMxRILs/geno_PLINK_Mx_chr6.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr7.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr7.txt" --header=2 --mode=2 --cross="Mo17" > BMxRILs/geno_PLINK_Mx_chr7.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr8.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr8.txt" --header=2 --mode=2 --cross="Mo17" > BMxRILs/geno_PLINK_Mx_chr8.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr9.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr9.txt" --header=2 --mode=2 --cross="Mo17" > BMxRILs/geno_PLINK_Mx_chr9.ped 
./myimpute_beta4.pl --dsnp="cleaned/clean_merged_chr10.dsnp" --bpfile="BMxRILs/tagging_Mx289_chr10.txt" --header=2 --mode=2 --cross="Mo17" > BMxRILs/geno_PLINK_Mx_chr10.ped 


###################################################################################
# NAM Diallels
###################################################################################
./myimpute_diallel.pl --diallel="Diallel/diallel_4impute.txt" --dsnp="cleaned/clean_merged_chr1.dsnp" --mode=2 --header=2 > Diallel/geno_PLINK_diallel_chr1.ped
./myimpute_diallel.pl --diallel="Diallel/diallel_4impute.txt" --dsnp="cleaned/clean_merged_chr2.dsnp" --mode=2 --header=2 > Diallel/geno_PLINK_diallel_chr2.ped
./myimpute_diallel.pl --diallel="Diallel/diallel_4impute.txt" --dsnp="cleaned/clean_merged_chr3.dsnp" --mode=2 --header=2 > Diallel/geno_PLINK_diallel_chr3.ped
./myimpute_diallel.pl --diallel="Diallel/diallel_4impute.txt" --dsnp="cleaned/clean_merged_chr4.dsnp" --mode=2 --header=2 > Diallel/geno_PLINK_diallel_chr4.ped
./myimpute_diallel.pl --diallel="Diallel/diallel_4impute.txt" --dsnp="cleaned/clean_merged_chr5.dsnp" --mode=2 --header=2 > Diallel/geno_PLINK_diallel_chr5.ped
./myimpute_diallel.pl --diallel="Diallel/diallel_4impute.txt" --dsnp="cleaned/clean_merged_chr6.dsnp" --mode=2 --header=2 > Diallel/geno_PLINK_diallel_chr6.ped
./myimpute_diallel.pl --diallel="Diallel/diallel_4impute.txt" --dsnp="cleaned/clean_merged_chr7.dsnp" --mode=2 --header=2 > Diallel/geno_PLINK_diallel_chr7.ped
./myimpute_diallel.pl --diallel="Diallel/diallel_4impute.txt" --dsnp="cleaned/clean_merged_chr8.dsnp" --mode=2 --header=2 > Diallel/geno_PLINK_diallel_chr8.ped
./myimpute_diallel.pl --diallel="Diallel/diallel_4impute.txt" --dsnp="cleaned/clean_merged_chr9.dsnp" --mode=2 --header=2 > Diallel/geno_PLINK_diallel_chr9.ped
./myimpute_diallel.pl --diallel="Diallel/diallel_4impute.txt" --dsnp="cleaned/clean_merged_chr10.dsnp" --mode=2 --header=2 > Diallel/geno_PLINK_diallel_chr10.ped


###### merge PLINK ----
cat IBM_impute/geno_PLINK_ibm136_chr1.PLINK >> geno_chr1.PLINK
cat BMxRILs/geno_PLINK_Bx_chr1.ped >> geno_chr1.PLINK
cat BMxRILs/geno_PLINK_Mx_chr1.ped >> geno_chr1.PLINK
cat Diallel/geno_PLINK_diallel_chr1.ped >> geno_chr1.PLINK

















###### merge GenSel
tail -n +2 IBM_impute/geno_IBM136_chr1.GenSel >> geno_meta_chr1.GenSel
tail -n +2 BMxRILs/geno_Bx692_chr1.GenSel >> geno_meta_chr1.GenSel
tail -n +2 BmxRILs/geno_Mx289_chr1.GenSel >> geno_meta_chr1.GenSel
tail -n +2 ~/Documents/Heterosis_GWAS/SNP/geno_diallel_GenSel_061312_chr1 >> geno_meta_chr1.GenSel

tail -n +2 IBM_impute/geno_IBM136_chr2.GenSel >> geno_meta_chr2.GenSel
tail -n +2 BMxRILs/geno_Bx692_chr2.GenSel >> geno_meta_chr2.GenSel
tail -n +2 BmxRILs/geno_Mx289_chr2.GenSel >> geno_meta_chr2.GenSel
tail -n +2 ~/Documents/Heterosis_GWAS/SNP/geno_diallel_GenSel_061312_chr2 >> geno_meta_chr2.GenSel

tail -n +2 IBM_impute/geno_IBM136_chr3.GenSel >> geno_meta_chr3.GenSel
tail -n +2 BMxRILs/geno_Bx692_chr3.GenSel >> geno_meta_chr3.GenSel
tail -n +2 BmxRILs/geno_Mx289_chr3.GenSel >> geno_meta_chr3.GenSel
tail -n +2 ~/Documents/Heterosis_GWAS/SNP/geno_diallel_GenSel_061312_chr3 >> geno_meta_chr3.GenSel

tail -n +2 IBM_impute/geno_IBM136_chr4.GenSel >> geno_meta_chr4.GenSel
tail -n +2 BMxRILs/geno_Bx692_chr4.GenSel >> geno_meta_chr4.GenSel
tail -n +2 BmxRILs/geno_Mx289_chr4.GenSel >> geno_meta_chr4.GenSel
tail -n +2 ~/Documents/Heterosis_GWAS/SNP/geno_diallel_GenSel_061312_chr4 >> geno_meta_chr4.GenSel

tail -n +2 IBM_impute/geno_IBM136_chr5.GenSel >> geno_meta_chr5.GenSel
tail -n +2 BMxRILs/geno_Bx692_chr5.GenSel >> geno_meta_chr5.GenSel
tail -n +2 BmxRILs/geno_Mx289_chr5.GenSel >> geno_meta_chr5.GenSel
tail -n +2 ~/Documents/Heterosis_GWAS/SNP/geno_diallel_GenSel_061312_chr5 >> geno_meta_chr5.GenSel

tail -n +2 IBM_impute/geno_IBM136_chr6.GenSel >> geno_meta_chr6.GenSel
tail -n +2 BMxRILs/geno_Bx692_chr6.GenSel >> geno_meta_chr6.GenSel
tail -n +2 BmxRILs/geno_Mx289_chr6.GenSel >> geno_meta_chr6.GenSel
tail -n +2 ~/Documents/Heterosis_GWAS/SNP/geno_diallel_GenSel_061312_chr6 >> geno_meta_chr6.GenSel

tail -n +2 IBM_impute/geno_IBM136_chr7.GenSel >> geno_meta_chr7.GenSel
tail -n +2 BMxRILs/geno_Bx692_chr7.GenSel >> geno_meta_chr7.GenSel
tail -n +2 BmxRILs/geno_Mx289_chr7.GenSel >> geno_meta_chr7.GenSel
tail -n +2 ~/Documents/Heterosis_GWAS/SNP/geno_diallel_GenSel_061312_chr7 >> geno_meta_chr7.GenSel

tail -n +2 IBM_impute/geno_IBM136_chr8.GenSel >> geno_meta_chr8.GenSel
tail -n +2 BMxRILs/geno_Bx692_chr8.GenSel >> geno_meta_chr8.GenSel
tail -n +2 BmxRILs/geno_Mx289_chr8.GenSel >> geno_meta_chr8.GenSel
tail -n +2 ~/Documents/Heterosis_GWAS/SNP/geno_diallel_GenSel_061312_chr8 >> geno_meta_chr8.GenSel

tail -n +2 IBM_impute/geno_IBM136_chr9.GenSel >> geno_meta_chr9.GenSel
tail -n +2 BMxRILs/geno_Bx692_chr9.GenSel >> geno_meta_chr9.GenSel
tail -n +2 BmxRILs/geno_Mx289_chr9.GenSel >> geno_meta_chr9.GenSel
tail -n +2 ~/Documents/Heterosis_GWAS/SNP/geno_diallel_GenSel_061312_chr9 >> geno_meta_chr9.GenSel

tail -n +2 IBM_impute/geno_IBM136_chr10.GenSel >> geno_meta_chr10.GenSel
tail -n +2 BMxRILs/geno_Bx692_chr10.GenSel >> geno_meta_chr10.GenSel
tail -n +2 BmxRILs/geno_Mx289_chr10.GenSel >> geno_meta_chr10.GenSel
tail -n +2 ~/Documents/Heterosis_GWAS/SNP/geno_diallel_GenSel_061312_chr10 >> geno_meta_chr10.GenSel



