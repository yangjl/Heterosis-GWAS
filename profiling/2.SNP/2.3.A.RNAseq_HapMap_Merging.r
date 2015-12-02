# Jinliang Yang
# Merging the two set of SNP
# updated: 6/08/2012

### HapMap2v2
plink --bfile zmHapMap2v2_chr1_s2 --merge-list mergelist.txt --make-bed --out zmHapMap2v2_all_s2
########################################
# HapMap1 merge with RNA-seq
########################################

plink --bfile NAM.AGPv2.SNPs.pooled --exclude HapMap1_RNAseq.missnp --make-bed  --out NAM.AGPv2.SNPs.pooled_exclude8k
plink --bfile NAM.AGPv2.SNPs.pooled_exclude8k --bmerge zmhapmap1_v2_060712.bed zmhapmap1_v2_060712.bim zmhapmap1_v2_060712.fam --merge-mode 7 --allow-no-sex  --out HapMap1_RNAseq

########################################
# HapMap1 merge with HapMap2
########################################
plink --bfile zmhapmap1_v2_060712 --exclude HapMap2_HapMap1_all.missnp --make-bed  --out zmHapMap1v2_exclude20k

plink --bfile zmHapMap1v2_exclude20k --chr 1 --bmerge zmHapMap2v2_chr1_s2.bed zmHapMap2v2_chr1_s2.bim zmHapMap2v2_chr1_s2.fam --merge-mode 7 --allow-no-sex  --out HapMap2_HapMap1_chr1
plink --bfile zmHapMap1v2_exclude20k --chr 2 --bmerge zmHapMap2v2_chr2_s2.bed zmHapMap2v2_chr2_s2.bim zmHapMap2v2_chr2_s2.fam --merge-mode 7 --allow-no-sex  --out HapMap2_HapMap1_chr2
plink --bfile zmHapMap1v2_exclude20k --chr 3 --bmerge zmHapMap2v2_chr3_s2.bed zmHapMap2v2_chr3_s2.bim zmHapMap2v2_chr3_s2.fam --merge-mode 7 --allow-no-sex  --out HapMap2_HapMap1_chr3
plink --bfile zmHapMap1v2_exclude20k --chr 4 --bmerge zmHapMap2v2_chr4_s2.bed zmHapMap2v2_chr4_s2.bim zmHapMap2v2_chr4_s2.fam --merge-mode 7 --allow-no-sex  --out HapMap2_HapMap1_chr4
plink --bfile zmHapMap1v2_exclude20k --chr 5 --bmerge zmHapMap2v2_chr5_s2.bed zmHapMap2v2_chr5_s2.bim zmHapMap2v2_chr5_s2.fam --merge-mode 7 --allow-no-sex  --out HapMap2_HapMap1_chr5
plink --bfile zmHapMap1v2_exclude20k --chr 6 --bmerge zmHapMap2v2_chr6_s2.bed zmHapMap2v2_chr6_s2.bim zmHapMap2v2_chr6_s2.fam --merge-mode 7 --allow-no-sex  --out HapMap2_HapMap1_chr6
plink --bfile zmHapMap1v2_exclude20k --chr 7 --bmerge zmHapMap2v2_chr7_s2.bed zmHapMap2v2_chr7_s2.bim zmHapMap2v2_chr7_s2.fam --merge-mode 7 --allow-no-sex  --out HapMap2_HapMap1_chr7
plink --bfile zmHapMap1v2_exclude20k --chr 8 --bmerge zmHapMap2v2_chr8_s2.bed zmHapMap2v2_chr8_s2.bim zmHapMap2v2_chr8_s2.fam --merge-mode 7 --allow-no-sex  --out HapMap2_HapMap1_chr8
plink --bfile zmHapMap1v2_exclude20k --chr 9 --bmerge zmHapMap2v2_chr9_s2.bed zmHapMap2v2_chr9_s2.bim zmHapMap2v2_chr9_s2.fam --merge-mode 7 --allow-no-sex  --out HapMap2_HapMap1_chr9
plink --bfile zmHapMap1v2_exclude20k --chr 10 --bmerge zmHapMap2v2_chr10_s2.bed zmHapMap2v2_chr10_s2.bim zmHapMap2v2_chr10_s2.fam --merge-mode 7 --allow-no-sex  --out HapMap2_HapMap1_chr10

########################################
# merge all three files
########################################
plink --bfile zmHapMap1v2_exclude20k --chr 1 --make-bed --allow-no-sex --out zmHapMap1v2_exclude20k_chr1
plink --bfile zmHapMap1v2_exclude20k --chr 2 --make-bed --allow-no-sex --out zmHapMap1v2_exclude20k_chr2
plink --bfile zmHapMap1v2_exclude20k --chr 3 --make-bed --allow-no-sex --out zmHapMap1v2_exclude20k_chr3
plink --bfile zmHapMap1v2_exclude20k --chr 4 --make-bed --allow-no-sex --out zmHapMap1v2_exclude20k_chr4
plink --bfile zmHapMap1v2_exclude20k --chr 5 --make-bed --allow-no-sex --out zmHapMap1v2_exclude20k_chr5
plink --bfile zmHapMap1v2_exclude20k --chr 6 --make-bed --allow-no-sex --out zmHapMap1v2_exclude20k_chr6
plink --bfile zmHapMap1v2_exclude20k --chr 7 --make-bed --allow-no-sex --out zmHapMap1v2_exclude20k_chr7
plink --bfile zmHapMap1v2_exclude20k --chr 8 --make-bed --allow-no-sex --out zmHapMap1v2_exclude20k_chr8
plink --bfile zmHapMap1v2_exclude20k --chr 9 --make-bed --allow-no-sex --out zmHapMap1v2_exclude20k_chr9
plink --bfile zmHapMap1v2_exclude20k --chr 10 --make-bed --allow-no-sex --out zmHapMap1v2_exclude20k_chr10


plink --bfile NAM.AGPv2.SNPs.pooled_exclude46k --chr 1 --merge-list mergelist_chr1.txt --make-bed --allow-no-sex  --out RNAseq_HapMap1_HapMap2_chr1
plink --bfile NAM.AGPv2.SNPs.pooled_exclude46k --chr 2 --merge-list mergelist_chr2.txt --make-bed --allow-no-sex  --out RNAseq_HapMap1_HapMap2_chr2
plink --bfile NAM.AGPv2.SNPs.pooled_exclude46k --chr 3 --merge-list mergelist_chr3.txt --make-bed --allow-no-sex  --out RNAseq_HapMap1_HapMap2_chr3
plink --bfile NAM.AGPv2.SNPs.pooled_exclude46k --chr 4 --merge-list mergelist_chr4.txt --make-bed --allow-no-sex  --out RNAseq_HapMap1_HapMap2_chr4
plink --bfile NAM.AGPv2.SNPs.pooled_exclude46k --chr 5 --merge-list mergelist_chr5.txt --make-bed --allow-no-sex  --out RNAseq_HapMap1_HapMap2_chr5
plink --bfile NAM.AGPv2.SNPs.pooled_exclude46k --chr 6 --merge-list mergelist_chr6.txt --make-bed --allow-no-sex  --out RNAseq_HapMap1_HapMap2_chr6
plink --bfile NAM.AGPv2.SNPs.pooled_exclude46k --chr 7 --merge-list mergelist_chr7.txt --make-bed --allow-no-sex  --out RNAseq_HapMap1_HapMap2_chr7
plink --bfile NAM.AGPv2.SNPs.pooled_exclude46k --chr 8 --merge-list mergelist_chr8.txt --make-bed --allow-no-sex  --out RNAseq_HapMap1_HapMap2_chr8
plink --bfile NAM.AGPv2.SNPs.pooled_exclude46k --chr 9 --merge-list mergelist_chr9.txt --make-bed --allow-no-sex  --out RNAseq_HapMap1_HapMap2_chr9
plink --bfile NAM.AGPv2.SNPs.pooled_exclude46k --chr 10 --merge-list mergelist_chr10.txt --make-bed --allow-no-sex  --out RNAseq_HapMap1_HapMap2_chr10

########################################
# Filtering MAF 0.1 and missing rate 0.4
########################################
plink --bfile RNAseq_HapMap1_HapMap2_chr1 --transpose --maf 0.1 --geno 0.4 --recode --out cleaned_SNP_chr1 # 0.6
plink --bfile RNAseq_HapMap1_HapMap2_chr2 --transpose --maf 0.1 --geno 0.4 --recode --out cleaned_SNP_chr2
plink --bfile RNAseq_HapMap1_HapMap2_chr3 --transpose --maf 0.1 --geno 0.4 --recode --out cleaned_SNP_chr3
plink --bfile RNAseq_HapMap1_HapMap2_chr4 --transpose --maf 0.1 --geno 0.4 --recode --out cleaned_SNP_chr4
plink --bfile RNAseq_HapMap1_HapMap2_chr5 --transpose --maf 0.1 --geno 0.4 --recode --out cleaned_SNP_chr5
plink --bfile RNAseq_HapMap1_HapMap2_chr6 --transpose --maf 0.1 --geno 0.4 --recode --out cleaned_SNP_chr6
plink --bfile RNAseq_HapMap1_HapMap2_chr7 --transpose --maf 0.1 --geno 0.4 --recode --out cleaned_SNP_chr7
plink --bfile RNAseq_HapMap1_HapMap2_chr8 --transpose --maf 0.1 --geno 0.4 --recode --out cleaned_SNP_chr8
plink --bfile RNAseq_HapMap1_HapMap2_chr9 --transpose --maf 0.1 --geno 0.4 --recode --out cleaned_SNP_chr9
plink --bfile RNAseq_HapMap1_HapMap2_chr10 --transpose --maf 0.1 --geno 0.4 --recode --out cleaned_SNP_chr10

########################################
# calculate the missing rate and MAF
########################################
plink --bfile RNAseq_HapMap1_HapMap2_chr1 --missing --out RNAseq_HapMap1_HapMap2_chr1
plink --bfile RNAseq_HapMap1_HapMap2_chr2 --missing --out RNAseq_HapMap1_HapMap2_chr2
plink --bfile RNAseq_HapMap1_HapMap2_chr3 --missing --out RNAseq_HapMap1_HapMap2_chr3
plink --bfile RNAseq_HapMap1_HapMap2_chr4 --missing --out RNAseq_HapMap1_HapMap2_chr4
plink --bfile RNAseq_HapMap1_HapMap2_chr5 --missing --out RNAseq_HapMap1_HapMap2_chr5
plink --bfile RNAseq_HapMap1_HapMap2_chr6 --missing --out RNAseq_HapMap1_HapMap2_chr6
plink --bfile RNAseq_HapMap1_HapMap2_chr7 --missing --out RNAseq_HapMap1_HapMap2_chr7
plink --bfile RNAseq_HapMap1_HapMap2_chr8 --missing --out RNAseq_HapMap1_HapMap2_chr8
plink --bfile RNAseq_HapMap1_HapMap2_chr9 --missing --out RNAseq_HapMap1_HapMap2_chr9
plink --bfile RNAseq_HapMap1_HapMap2_chr10 --missing --out RNAseq_HapMap1_HapMap2_chr10

plink --bfile RNAseq_HapMap1_HapMap2_chr1 --freq --out RNAseq_HapMap1_HapMap2_chr1
plink --bfile RNAseq_HapMap1_HapMap2_chr2 --freq --out RNAseq_HapMap1_HapMap2_chr2
plink --bfile RNAseq_HapMap1_HapMap2_chr3 --freq --out RNAseq_HapMap1_HapMap2_chr3
plink --bfile RNAseq_HapMap1_HapMap2_chr4 --freq --out RNAseq_HapMap1_HapMap2_chr4
plink --bfile RNAseq_HapMap1_HapMap2_chr5 --freq --out RNAseq_HapMap1_HapMap2_chr5
plink --bfile RNAseq_HapMap1_HapMap2_chr6 --freq --out RNAseq_HapMap1_HapMap2_chr6
plink --bfile RNAseq_HapMap1_HapMap2_chr7 --freq --out RNAseq_HapMap1_HapMap2_chr7
plink --bfile RNAseq_HapMap1_HapMap2_chr8 --freq --out RNAseq_HapMap1_HapMap2_chr8
plink --bfile RNAseq_HapMap1_HapMap2_chr9 --freq --out RNAseq_HapMap1_HapMap2_chr9
plink --bfile RNAseq_HapMap1_HapMap2_chr10 --freq --out RNAseq_HapMap1_HapMap2_chr10

###--- Plot the results
chr1 <- read.table("RNAseq_HapMap1_HapMap2_chr1.lmiss", header=T)
chr2 <- read.table("RNAseq_HapMap1_HapMap2_chr2.lmiss", header=T)
chr3 <- read.table("RNAseq_HapMap1_HapMap2_chr3.lmiss", header=T)
chr4 <- read.table("RNAseq_HapMap1_HapMap2_chr4.lmiss", header=T)
chr5 <- read.table("RNAseq_HapMap1_HapMap2_chr5.lmiss", header=T)
chr6 <- read.table("RNAseq_HapMap1_HapMap2_chr6.lmiss", header=T)
chr7 <- read.table("RNAseq_HapMap1_HapMap2_chr7.lmiss", header=T)
chr8 <- read.table("RNAseq_HapMap1_HapMap2_chr8.lmiss", header=T)
chr9 <- read.table("RNAseq_HapMap1_HapMap2_chr9.lmiss", header=T)
chr10 <- read.table("RNAseq_HapMap1_HapMap2_chr10.lmiss", header=T)

lmiss <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10)

chr1 <- read.table("RNAseq_HapMap1_HapMap2_chr1.frq", header=T)
chr2 <- read.table("RNAseq_HapMap1_HapMap2_chr2.frq", header=T)
chr3 <- read.table("RNAseq_HapMap1_HapMap2_chr3.frq", header=T)
chr4 <- read.table("RNAseq_HapMap1_HapMap2_chr4.frq", header=T)
chr5 <- read.table("RNAseq_HapMap1_HapMap2_chr5.frq", header=T)
chr6 <- read.table("RNAseq_HapMap1_HapMap2_chr6.frq", header=T)
chr7 <- read.table("RNAseq_HapMap1_HapMap2_chr7.frq", header=T)
chr8 <- read.table("RNAseq_HapMap1_HapMap2_chr8.frq", header=T)
chr9 <- read.table("RNAseq_HapMap1_HapMap2_chr9.frq", header=T)
chr10 <- read.table("RNAseq_HapMap1_HapMap2_chr10.frq", header=T)

frq <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10)

par(mfrow=c(1,2))
hist(lmiss$F_MISS, main="Histogram of loci Missing Rate", xlab="missing rate", breaks=25, col="antiquewhite3");
hist(frq$MAF, main="Histogram of MAF", xlab="MAF", breaks=30, col="antiquewhite3")

subsnp1 <- subset(frq, MAF > 0.1)
subsnp2 <- subset(lmiss, F_MISS > 0.6)
subsnp <- merge(subsnp1, subsnp2, by="SNP")
dim(subsnp)
#[1] 1736352      10

write.table(subsnp, "pav_snp_1.7M.txt", sep="\t", quote=FALSE, row.names=FALSE)
 




