#Jinliang yang
#Purpose: Use python code to compute snp frq and missingness
#updated: 7.5.2014


setwd("~/DBcenter/VariationDB/merged/")

python ~/Documents/PyCodes/snpfrq_v1.0.py -i snp_merged_chr1 -s 5 -e 31 -o snp_merged_chr1_frq
python ~/Documents/PyCodes/snpfrq_v1.0.py -i snp_merged_chr2 -s 5 -e 31 -o snp_merged_chr2_frq
python ~/Documents/PyCodes/snpfrq_v1.0.py -i snp_merged_chr3 -s 5 -e 31 -o snp_merged_chr3_frq
python ~/Documents/PyCodes/snpfrq_v1.0.py -i snp_merged_chr4 -s 5 -e 31 -o snp_merged_chr4_frq
python ~/Documents/PyCodes/snpfrq_v1.0.py -i snp_merged_chr5 -s 5 -e 31 -o snp_merged_chr5_frq
python ~/Documents/PyCodes/snpfrq_v1.0.py -i snp_merged_chr6 -s 5 -e 31 -o snp_merged_chr6_frq
python ~/Documents/PyCodes/snpfrq_v1.0.py -i snp_merged_chr7 -s 5 -e 31 -o snp_merged_chr7_frq
python ~/Documents/PyCodes/snpfrq_v1.0.py -i snp_merged_chr8 -s 5 -e 31 -o snp_merged_chr8_frq
python ~/Documents/PyCodes/snpfrq_v1.0.py -i snp_merged_chr9 -s 5 -e 31 -o snp_merged_chr9_frq
python ~/Documents/PyCodes/snpfrq_v1.0.py -i snp_merged_chr10 -s 5 -e 31 -o snp_merged_chr10_frq


##test the codes
setwd("~/DBcenter/VariationDB/merged/")

tab5rows <- read.table("snp_merged_chr1", nrow=5, header=FALSE)
classes <- sapply(tab5rows, class)
#classes[2] <- "character"
dsf <- read.table("snp_merged_chr1", header = FALSE, colClasses = classes)

snp <- read.table("snp_merged_chr1", header=TRUE)
frq <- read.table("snp_merged_chr1_frq", header=TRUE)



hmp1 <- read.table("~/DBcenter/VariationDB/HapMap1/hapmap1V2_070214_chr1.dsf", header=T)
hmp2 <- read.table("~/DBcenter/VariationDB/HapMap2/maizeHapMap2V2_chr1.dsf", header=T)

subset(hmp1, snpid==snpid)
subset(hmp2, snpid==snpid)
subset(rnaseq, snpid==snpid)


