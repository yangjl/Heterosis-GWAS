## Jinliang Yang
## July 17th, 2014
## plot the phenotypic traits


readPheno <- function(infile="akw_v2_GenSel_fullset.txt"){
  setwd("~/Documents/Heterosis_GWAS/HGWAS_proj/data/")
  pheno <- read.table(infile, header=TRUE)
  names(pheno)[3:5] <- c("trait", "pop", "FID");
  pheno$trait <- names(pheno)[2]
  names(pheno)[2] <- "value"
  message(sprintf("[ %s ] genotypes", nrow(pheno)))
  return(pheno)
}


######
krn <- readPheno(infile="krn_v2_GenSel_fullset.txt")
akw <- readPheno(infile="akw_v2_GenSel_fullset.txt")
cd <- readPheno(infile="cd_v2_GenSel_fullset.txt")
cl <- readPheno(infile="cl_v2_GenSel_fullset.txt")
cw <- readPheno(infile="cw_v2_GenSel_fullset.txt")
kc <- readPheno(infile="kc_v2_GenSel_fullset.txt")
tkw <- readPheno(infile="tkw_v2_GenSel_fullset.txt")

trait <- rbind(krn, akw, cd, cl, cw, kc, tkw)
write.table(trait, "~/Documents/Heterosis_GWAS/HGWAS_proj/reports/S1_pheno.txt",
            sep="\t", row.names=FALSE, quote=FALSE)

#######################
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
write.table(pmat, "~/Documents/Heterosis_GWAS/HGWAS_proj/cache/pheno_matrix.txt",
            sep="\t", row.names=FALSE, quote=FALSE)
#################################################

pmat <- read.table("~/Documents/Heterosis_GWAS/HGWAS_proj/cache/pheno_matrix.txt", header=T)


plotKRN <- function(){
  krnall$norv <- scale(krnall$value)
  plot(density(subset(krnall, pop=="NAMRIL")$norv), col="black", lwd=4, ylim=c(0, 0.8), bty="n",
       main="Kernel Row Number", xlab="KRN")
  lines(density(subset(krnall, pop=="Diallel")$norv), col="gold", lwd=4)
  lines(density(subset(krnall, pop=="BxRIL")$norv), col="blue", lwd=4)
  lines(density(subset(krnall, pop=="MxRIL")$norv), col="red", lwd=4)
  #abline(v=17.1, col="blue", lty=2, lwd=2)
  temp <- legend("topright", legend =c(" ", " ", " ", " "), title="Population", 
                 text.width=strwidth("MxRIL"),
                 lty=1, lwd=3, col=c("black", "gold", "blue", "red"), 
                 xjust=1, yjust=1)
  text(temp$rect$left + temp$rect$w, temp$text$y,
       c("RIL", "Diallel", "BxRIL", "MxRIL"), pos=2)  
}

###########
denplot <- function(tcol="KRN", ...){
  par(mar=c(2,2,3,2))
  myp <- subset(trait, trait==tcol)
  myp$norv <- scale(myp$value)
  plot(density(subset(myp, pop=="NAMRIL")$norv), col="black", lwd=4, bty="n",
       main=tcol, xlab="", ...)
  lines(density(subset(myp, pop=="Diallel")$norv), col="gold", lwd=4)
  lines(density(subset(myp, pop=="BxRIL")$norv), col="blue", lwd=4)
  lines(density(subset(myp, pop=="MxRIL")$norv), col="red", lwd=4)
}

#######
#nf <- layout(matrix(c(1,1,2,3,4,1,1,5,6,7), 2, 5, byrow = TRUE))
#layout.show(nf)



pdf("~/Documents/Heterosis_GWAS/HGWAS_proj/reports/S.F1_pheno.pdf", width=10, height=5)
layout(matrix(c(1,1,2,3,4,1,1,5,6,7), 2, 5, byrow = TRUE))
#par(mar=c(0,0,0,0))
plotKRN()
denplot(tcol="CD", xlim=c(-4, 4), ylim=c(0, 1))
denplot(tcol="AKW", xlim=c(-4, 5), ylim=c(0, 1))
denplot(tcol="CL", xlim=c(-4, 5), ylim=c(0, 0.8))
denplot(tcol="CW", xlim=c(-4, 8), ylim=c(0, 0.8))
denplot(tcol="KC", xlim=c(-4, 6), ylim=c(0, 0.8))
denplot(tcol="TKW", xlim=c(-4, 6), ylim=c(0, 1.2))
dev.off()


