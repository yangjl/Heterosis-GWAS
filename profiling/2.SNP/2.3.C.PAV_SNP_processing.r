# Jinliang Yang
# Merging the two set of SNP
# updated: 6/08/2012

setwd("/Users/yangjl/Documents/GWAS2_KRN/SNP/merged/all_merged")

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

write.table(subsnp[,1], "pav_snp_1.7M.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
 
########################################
plink --bfile RNAseq_HapMap1_HapMap2_chr1 --extract pav_snp_1.7M.txt --make-bed --out pav_SNP_chr1
plink --bfile RNAseq_HapMap1_HapMap2_chr2 --extract pav_snp_1.7M.txt --make-bed --out pav_SNP_chr2
plink --bfile RNAseq_HapMap1_HapMap2_chr3 --extract pav_snp_1.7M.txt --make-bed --out pav_SNP_chr3
plink --bfile RNAseq_HapMap1_HapMap2_chr4 --extract pav_snp_1.7M.txt --make-bed --out pav_SNP_chr4
plink --bfile RNAseq_HapMap1_HapMap2_chr5 --extract pav_snp_1.7M.txt --make-bed --out pav_SNP_chr5
plink --bfile RNAseq_HapMap1_HapMap2_chr6 --extract pav_snp_1.7M.txt --make-bed --out pav_SNP_chr6
plink --bfile RNAseq_HapMap1_HapMap2_chr7 --extract pav_snp_1.7M.txt --make-bed --out pav_SNP_chr7
plink --bfile RNAseq_HapMap1_HapMap2_chr8 --extract pav_snp_1.7M.txt --make-bed --out pav_SNP_chr8
plink --bfile RNAseq_HapMap1_HapMap2_chr9 --extract pav_snp_1.7M.txt --make-bed --out pav_SNP_chr9
plink --bfile RNAseq_HapMap1_HapMap2_chr10 --extract pav_snp_1.7M.txt --make-bed --out pav_SNP_chr10

##### merge
plink --bfile pav_SNP_chr1 --merge-list mergelist_pav_snp.txt --allow-no-sex --transpose --recode --out pav_SNP_all
#### R
pavsnp <- read.table("pav_SNP_all.tped", header=FALSE)
pav1 <- pavsnp[,5:58]
pav1[pav1!=0] <- "A"; #present
pav1[pav1==0] <- "T"; #absent
pavsnp[,5:58] <- pav1;
write.table(pavsnp, "pav_SNP_all.tped", sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE);
###############
plink --tfile pav_SNP_all --maf 0.1 --recode --transpose --out pav_SNP_all
#1435003 SNPs

tped2dsnp <- function(tpedfile="cleaned_SNP_chr1.tped", outfile="clean_merged_chr1.dsnp"){
	
	nm <- c("chr", "rs","allele", "pos","B73", "Z001","Z002", "Z003","Z004", 
			"Z005","Z006","Z007","Z008","Z009","Z010",
			"Z011","Z012","Z013","Z014","Z015","Z016",
			"Z017","Z018","Z019","Z020","Z021","Z022",
			"Z023","Z024","Z025","Z026")
	
	tped1 <- read.table(tpedfile, header=FALSE)
	tped1 <- tped1[, c(1:4, 5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57)]
	names(tped1) <- nm
	tem <- tped1[, 5:31]
	tem <- apply(tem, 2, as.character)
	tem[tem==0] <- "N"
	tped1[,5:31] <- tem
	write.table(tped1, outfile, sep="\t", row.names=FALSE, quote=FALSE)
	return(dim(tped1))
}

tped2dsnp(tpedfile="pav_SNP_all.tped", outfile="pav_SNP_all.dsnp")










