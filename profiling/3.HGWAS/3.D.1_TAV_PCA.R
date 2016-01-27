### Jinliang Yang
### conduct PCA 

tavsnp <- read.csv("largedata/SNP/TAV_recoded.csv")

mx <- as.matrix(tavsnp[, -1])
sampleid <- read.table("largedata/SNP/TKW_chr9-covars.sample", header=TRUE)


snpgdsCreateGeno("largedata/lcache/tavsnp.gds", genmat=mx,
                 sample.id= sampleid$ID_1[-1], snp.id= as.character(tavsnp$V1), 
                 snp.chromosome= gsub("_.*", "", tavsnp$V1),
                 snp.position= gsub(".*_", "", tavsnp$V1), snp.allele=NULL, snpfirstdim=TRUE)

library(SNPRelate)
(genofile <- snpgdsOpen("largedata/lcache/tavsnp.gds"))

pca <- snpgdsPCA(genofile, num.thread=2)
snpgdsClose(genofile)
# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
tab$sample.id <- as.character(tab$sample.id)
tab$sub1 <- gsub("_.*", "", tab$sample.id)
tab$sub2 <- gsub(".*_", "", tab$sample.id)
tab$pop <- "NAM"
idx1 <- grep("^B73x", tab$sub1)
idx2 <- grep("^Mo17x", tab$sub1)
idx3 <- grep("^1$", tab$sub2)

tab[idx3, ]$pop <- "diallel"
tab[idx1[!idx1%in%idx3], ]$pop <- "B73"
tab[idx2, ]$pop <- "Mo17"

z9 <- grep("Z009", tab$sample.id)

pdf("graphs/Test_PCA.pdf")
plot(tab[-z9,]$EV2, tab[-z9,]$EV1, col=as.integer(as.factor(tab[-z9,]$pop)), xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomleft", legend=levels(as.factor(tab$pop)), pch="o", col=1:nlevels(as.factor(tab$pop)))
dev.off()

pc.percent <- pca$varprop*100
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=as.integer(as.factor(tab$pop)), labels=lbls)

