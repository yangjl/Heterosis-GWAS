### Jinliang Yang
### 9/4/2014
### get dsf5 for diallel dominant imputationob2


snp14m <- read.table("/home/NSF-SAM-GS/vardb/NAM-diallel/chrall_14m.txt", header=TRUE)
##### SNPs
ob2 <- load("~/Documents/Heterosis_GWAS/HGWAS_proj/cache/bayes.RData")

head(snp14m)
#################
write_dsf <- function(snp14m, snp2=krn2,
                      outfile="/home/NSF-SAM-GS/HGWAS/SNPdiallel/snpkrn_1m.dsf5"){
  
  mysnp <- subset(snp14m, snpid %in% snp2$ID) #[1] 1153214      35
  message(sprintf("[ %s ] SNPs!", nrow(mysnp)))
  write.table(mysnp, outfile, row.names=FALSE,
              quote=FALSE, sep="\t")  
}
#########
write_dsf(snp14m, snp2=krn2, outfile="/home/NSF-SAM-GS/HGWAS/SNPdiallel/snpkrn_1m.dsf5")
write_dsf(snp14m, snp2=cd2, outfile="/home/NSF-SAM-GS/HGWAS/SNPdiallel/snpcd_1m.dsf5")
write_dsf(snp14m, snp2=cw2, outfile="/home/NSF-SAM-GS/HGWAS/SNPdiallel/snpcw_1m.dsf5")
write_dsf(snp14m, snp2=cl2, outfile="/home/NSF-SAM-GS/HGWAS/SNPdiallel/snpcl_1m.dsf5")
write_dsf(snp14m, snp2=akw2, outfile="/home/NSF-SAM-GS/HGWAS/SNPdiallel/snpakw_1m.dsf5")
write_dsf(snp14m, snp2=tkw2, outfile="/home/NSF-SAM-GS/HGWAS/SNPdiallel/snptkw_1m.dsf5")
write_dsf(snp14m, snp2=kc2, outfile="/home/NSF-SAM-GS/HGWAS/SNPdiallel/snpkc_1m.dsf5")

######## pedigree
setwd("~/Documents/Heterosis_GWAS/HGWAS_proj/data/")
pheno <- read.table("krn_v2_GenSel_fullset.txt", header=TRUE)
ped <- subset(pheno, pheno[,4] == "Diallel")
ped$p1 <- gsub("x.*", "", ped$Genotype)
ped$p2 <- gsub(".*x", "", ped$Genotype)
ped <- ped[, c("p1", "p2", "Genotype")]
names(ped) <- c("p1", "p2", "f1")
write.table(ped, "/home/NSF-SAM-GS/HGWAS/SNPdiallel/ped221.txt", sep="\t",
            row.names=FALSE, quote=FALSE)


### dominant
python impute4diallel_v1.2.py -h
cd /home/NSF-SAM-GS/HGWAS/SNPdiallel/

python ~/Documents/PyCodes/impute4diallel/impute4diallel_v1.2.py \
-d ped221.txt -i snpkrn_1m.dsf5 -o snpkrn_1m_dom_gensel -s 9 -e 35 --header yes -m 4
python ~/Documents/PyCodes/impute4diallel/impute4diallel_v1.2.py \
-d ped221.txt -i snpcd_1m.dsf5 -o snpcd_1m_dom_gensel -s 9 -e 35 --header yes -m 4
python ~/Documents/PyCodes/impute4diallel/impute4diallel_v1.2.py \
-d ped221.txt -i snpcw_1m.dsf5 -o snpcw_1m_dom_gensel -s 9 -e 35 --header yes -m 4
python ~/Documents/PyCodes/impute4diallel/impute4diallel_v1.2.py \
-d ped221.txt -i snpcl_1m.dsf5 -o snpcl_1m_dom_gensel -s 9 -e 35 --header yes -m 4
python ~/Documents/PyCodes/impute4diallel/impute4diallel_v1.2.py \
-d ped221.txt -i snpakw_1m.dsf5 -o snpakw_1m_dom_gensel -s 9 -e 35 --header yes -m 4
python ~/Documents/PyCodes/impute4diallel/impute4diallel_v1.2.py \
-d ped221.txt -i snptkw_1m.dsf5 -o snptkw_1m_dom_gensel -s 9 -e 35 --header yes -m 4
python ~/Documents/PyCodes/impute4diallel/impute4diallel_v1.2.py \
-d ped221.txt -i snpkc_1m.dsf5 -o snpkc_1m_dom_gensel -s 9 -e 35 --header yes -m 4



##############
python ~/Documents/PyCodes/impute4diallel/impute4diallel_v1.2.py \
-d ped221.txt -i snpkrn_1m.dsf5 -o snpkrn_1m_add_gensel -s 9 -e 35 --header yes -m 0
python ~/Documents/PyCodes/impute4diallel/impute4diallel_v1.2.py \
-d ped221.txt -i snpcd_1m.dsf5 -o snpcd_1m_add_gensel -s 9 -e 35 --header yes -m 0
python ~/Documents/PyCodes/impute4diallel/impute4diallel_v1.2.py \
-d ped221.txt -i snpcw_1m.dsf5 -o snpcw_1m_add_gensel -s 9 -e 35 --header yes -m 0
python ~/Documents/PyCodes/impute4diallel/impute4diallel_v1.2.py \
-d ped221.txt -i snpcl_1m.dsf5 -o snpcl_1m_add_gensel -s 9 -e 35 --header yes -m 0
python ~/Documents/PyCodes/impute4diallel/impute4diallel_v1.2.py \
-d ped221.txt -i snpakw_1m.dsf5 -o snpakw_1m_add_gensel -s 9 -e 35 --header yes -m 0
python ~/Documents/PyCodes/impute4diallel/impute4diallel_v1.2.py \
-d ped221.txt -i snptkw_1m.dsf5 -o snptkw_1m_add_gensel -s 9 -e 35 --header yes -m 0
python ~/Documents/PyCodes/impute4diallel/impute4diallel_v1.2.py \
-d ped221.txt -i snpkc_1m.dsf5 -o snpkc_1m_add_gensel -s 9 -e 35 --header yes -m 0










