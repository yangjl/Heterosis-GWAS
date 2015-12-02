# Jinliang Yang
# purpose: Merged SNP Imputation
# Last update: 7.17.2012


# GenSel
###################################################################################
# IBM RILs
###################################################################################
./myimpute_beta4.pl --dsnp="all_merged/pav_SNP_all.dsnp" --bpfile="IBM_impute/tagging_IBM136_all.txt" --header=1 --mode=1 --cross="RIL" > IBM_impute/pavsnp_ibm136_GenSel

###################################################################################
# NAM RILs
###################################################################################
./myimpute_beta4.pl --dsnp="all_merged/pav_SNP_all.dsnp" --bpfile="cleaned/tagging_gbs50k_all.txt" --header=2 --mode=1 --cross="RIL" > pavsnp_NAMRIL_GenSel

###################################################################################
# B73 x NAM RILs
###################################################################################
./myimpute_beta4.pl --dsnp="all_merged/pav_SNP_all.dsnp" --bpfile="BMxRILs/tagging_Bx692_all.txt" --header=2 --mode=1 --cross="B73" > BMxRILs/pavsnp_Bx_GenSel

###################################################################################
# Mo17 x NAM RILs
###################################################################################
./myimpute_beta4.pl --dsnp="all_merged/pav_SNP_all.dsnp" --bpfile="BMxRILs/tagging_Mx289_all.txt" --header=2 --mode=1 --cross="Mo17" > BMxRILs/pavsnp_Mx_GenSel 

###################################################################################
# NAM Diallels
###################################################################################
./myimpute_diallel.pl --diallel="Diallel/diallel_4impute.txt" --dsnp="all_merged/pav_SNP_all.dsnp" --header=2 --mode=1 > Diallel/pavsnp_diallel_GenSel


###### merge PLINK ----
mv IBM_impute/pavsnp_ibm136_GenSel geno_pavsnp_GenSel
cat pavsnp_NAMRIL_GenSel >> geno_pavsnp_GenSel
cat BMxRILs/pavsnp_Bx_GenSel >> geno_pavsnp_GenSel
cat BMxRILs/pavsnp_Mx_GenSel >> geno_pavsnp_GenSel
cat Diallel/pavsnp_diallel_GenSel >> geno_pavsnp_GenSel

./newgen2bin -i geno_pavsnp_GenSel





