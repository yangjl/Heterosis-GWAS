### Jinliang Yang
### January 5th, 2016

#### phenotypic data
sp <- read.table("largedata/SNP/TKW_chr9-covars.sample", header=TRUE)
sp <- sp[-1, 1:6]
names(sp)[1] <- "sampleid"
sp$order <- 1:nrow(sp)
naidx <- subset(sp, is.na(pop))$order

idx1 <- naidx[grep("^Z", sp[naidx, ]$sampleid)]
idx2 <- naidx[grep("^B73", sp[naidx, ]$sampleid)]

sp[sp$order %in% idx1, ]$pop <- 1
sp[sp$order %in% idx1, ]$FID <- c(5,6,10, 11,17,18,26)
sp[sp$order %in% idx2, ]$pop <- 2
sp[sp$order %in% idx2, ]$FID <- 17

sp$sampleid <- gsub("_", "E", sp$sampleid )
sp$sampleid  <- gsub("E1$", "", sp$sampleid )

pheno <- read.table("cache/pheno_matrix.txt", header=TRUE)
sp <- merge(sp[, -4], pheno[, -2:-3], by.x="sampleid", by.y="Genotype", all.x=TRUE)
sp <- sp[order(sp$order),]

for(i in 7:13){
    sp[is.na(sp[, i]), i] <- -999
}

write.table(sp, "largedata/pheno/pheno_sample.txt", sep="\t", row.names=FALSE, quote=FALSE)



### genotypic data
geno <- read.csv("largedata/SNP/TAV_recoded.csv")
names(geno) <- c("sampleid", sp$sampleid)
geno$sampleid <- paste0("S", geno$sampleid)
geno <- t(geno)
write.table(geno, "largedata/SNP/geno_gblup.txt", sep="\t", row.names=T, col.names=F, quote=FALSE)

### map file
geno <- read.csv("largedata/SNP/TAV_recoded.csv")

map <- geno[, 1:3]
names(map) <- c("snpid", "chr", "pos")
map$chr <- gsub("_.*", "", map$snpid)
map$pos <- gsub(".*_", "", map$snpid)
write.table(geno, "largedata/SNP/geno_gblup_map.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)

