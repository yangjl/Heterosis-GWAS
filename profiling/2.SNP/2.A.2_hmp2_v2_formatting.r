#Jinliang yang
#Purpose: formatting the HapMap2 V2 SNPs for imputation
#start: 6.6.2012
# updated: 7/1/2014

#Note: Running on server 129.186.85.7

# Formatting Function:###################################################
# 1. [remove the B73 != Ref SNP] do not do this anymore!
# 2. remove the duplicated SNP
# 3. [change the coding of -/+ to A/T];
# 4. change the IUPAC coding to "N"
# 5. Output to DSF (density snp format)
######################################################################

source("~/Documents/Heterosis_GWAS/HGWAS_proj/lib/chr2ped.R")
setwd("~/DBcenter/VariationDB/HapMap2/")

####
chr2dsf(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr2.hmp.txt", output="maizeHapMap2V2_chr2.dsf")

chr2dsf(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr1.hmp.txt", output="maizeHapMap2V2_chr1.dsf")
chr2dsf(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr3.hmp.txt", output="maizeHapMap2V2_chr3.dsf")
chr2dsf(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr4.hmp.txt", output="maizeHapMap2V2_chr4.dsf")
chr2dsf(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr5.hmp.txt", output="maizeHapMap2V2_chr5.dsf")
chr2dsf(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr6.hmp.txt", output="maizeHapMap2V2_chr6.dsf")
chr2dsf(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr7.hmp.txt", output="maizeHapMap2V2_chr7.dsf")
chr2dsf(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr8.hmp.txt", output="maizeHapMap2V2_chr8.dsf")
chr2dsf(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr9.hmp.txt", output="maizeHapMap2V2_chr9.dsf")
chr2dsf(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr10.hmp.txt", output="maizeHapMap2V2_chr10.dsf")


