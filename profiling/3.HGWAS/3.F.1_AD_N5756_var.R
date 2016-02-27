### Jinliang Yang
### Feb 26th, 2016

sp <- read.table("largedata/pheno/pheno_sample.txt", header=TRUE)
sp$sampleid <- as.character(sp$sampleid)
################################################################################################
### genotypic data
geno <- read.csv("largedata/SNP/TAV_recoded_N5756.csv")
names(geno) <- c("sampleid", sp$sampleid)
geno$sampleid <- paste0("S", geno$sampleid)
tgeno <- t(geno)
write.table(tgeno, "largedata/SNP/geno_gblup_new.txt", sep="\t", row.names=T, col.names=F, quote=FALSE)

### map file
geno <- read.csv("largedata/SNP/TAV_recoded_N5756.csv")
map <- geno[, 1:3]
names(map) <- c("snpid", "chr", "pos")
map$chr <- gsub("_.*", "", map$snpid)
map$pos <- gsub(".*_", "", map$snpid)
write.table(map, "largedata/SNP/geno_gblup_map_new.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)

