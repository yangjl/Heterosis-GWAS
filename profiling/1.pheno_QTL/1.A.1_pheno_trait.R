## Jinliang Yang
## July 17th, 2014
## plot the phenotypic traits
# updated: 12.2.2015


readPheno <- function(infile="akw_v2_GenSel_fullset.txt"){
  pheno <- read.table(infile, header=TRUE)
  names(pheno)[3:5] <- c("trait", "pop", "FID");
  pheno$trait <- names(pheno)[2]
  names(pheno)[2] <- "value"
  message(sprintf("[ %s ] genotypes", nrow(pheno)))
  return(pheno)
}


######
krn <- readPheno(infile="data/krn_v2_GenSel_fullset.txt")
akw <- readPheno(infile="data/akw_v2_GenSel_fullset.txt")
cd <- readPheno(infile="data/cd_v2_GenSel_fullset.txt")
cl <- readPheno(infile="data/cl_v2_GenSel_fullset.txt")
cw <- readPheno(infile="data/cw_v2_GenSel_fullset.txt")
kc <- readPheno(infile="data/kc_v2_GenSel_fullset.txt")
tkw <- readPheno(infile="data/tkw_v2_GenSel_fullset.txt")

trait <- rbind(krn, akw, cd, cl, cw, kc, tkw)
write.table(trait, "reports/S1_pheno.txt",
            sep="\t", row.names=FALSE, quote=FALSE)

#######################
# cache the phenotype matrix
names(krn)[2] <- "KRN"
names(cd)[2] <- "CD"
names(akw)[2] <- "AKW"
names(cl)[2] <- "CL"
names(cw)[2] <- "CW"
names(kc)[2] <- "KC"
names(tkw)[2] <- "TKW"
pmat <- merge(cd[, c(1, 4:5, 2)], krn[,1:2], by="Genotype", all=TRUE)
pmat <- merge(pmat, akw[, 1:2], by="Genotype", all=TRUE)
pmat <- merge(pmat, cl[, 1:2], by="Genotype", all=TRUE)
pmat <- merge(pmat, cw[, 1:2], by="Genotype", all=TRUE)
pmat <- merge(pmat, kc[, 1:2], by="Genotype", all=TRUE)
pmat <- merge(pmat, tkw[, 1:2], by="Genotype", all=TRUE)
pmat <- pmat[order(pmat$pop, pmat$Genotype),]
write.table(pmat, "cache/pheno_matrix.txt",sep="\t", row.names=FALSE, quote=FALSE)


