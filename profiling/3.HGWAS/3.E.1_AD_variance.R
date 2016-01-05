### Jinliang Yang
### January 5th, 2016

#### phenotypic data
sp <- read.table("largedata/SNP/TKW_chr9-covars.sample", header=TRUE)
sp <- sp[-1, 1:6]
names(sp)[1] <- "sampleid"

sp$sampleid <- gsub("_", "E", sp$sampleid )
sp$order <- 1:nrow(sp)
sp$sampleid  <- gsub("E1$", "", sp$sampleid )

pheno <- read.table("cache/pheno_matrix.txt", header=TRUE)
sp <- merge(sp[, -4], pheno[, -2:-3], by.x="ID_1", by.y="Genotype", all.x=TRUE)
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


subset(sp[, -2:-3], pop==2)