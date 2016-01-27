### Jinliang Yang
### GenSel
### 8/13/2014

readGenSel <- function(infile="CD_run41000.mrkRes1"){

  tab5rows <- read.table(infile, nrow=5, header=TRUE)
  classes <- sapply(tab5rows, class)
  
  pval <- read.table(infile, header=TRUE, colClasses=classes)
  pval <- pval[, 1:6]
  pval$chr <- gsub("_.*", "", pval$ID)
  pval$pos <- gsub(".*_", "", pval$ID)
  
  message(sprintf("TOTAL SNPs [ %s ]", nrow(pval)))
  return(pval)
}

########
setwd("~/Documents/Heterosis_GWAS/Method/GenSel")
tkw2 <- readGenSel(infile="TKW_run41000.mrkRes1")
krn2 <- readGenSel(infile="KRN_run41000.mrkRes1")
akw2 <- readGenSel(infile="AKW_run41000.mrkRes1")
kc2 <- readGenSel(infile="KC_run41000.mrkRes1")
cd2 <- readGenSel(infile="CD_run41000.mrkRes1")
cl2 <- readGenSel(infile="CL_run41000.mrkRes1")
cw2 <- readGenSel(infile="CW_run41000.mrkRes1")

source("~/Documents/Rcodes/save.append.R")
save.append(list=c("tkw2", "akw2", "kc2", "cl2", "cd2", "cw2", "krn2"), 
            file="~/Documents/Heterosis_GWAS/HGWAS_proj/cache/bayes.RData",
            description="")


#### quick check the SNPTEST results
source("~/Documents/Rcodes/quickMHTplot.R")

ob2 <- load("~/Documents/Heterosis_GWAS/HGWAS_proj/cache/bayes.RData")
qplot2 <- function(pval=krn2, ...){
  
  pval <- subset(pval, ModelFreq > 0.0001)
  pval$chr <- as.numeric(as.character(pval$chr))
  pval$pos <- as.numeric(as.character(pval$pos))
  message(sprintf("# of ModelFreq > 0.001 [ %s ]", nrow(pval)))
  quickMHTplot(res=pval, cex=.3, pch=19, col=rep(c("slateblue", "cyan4"), 5), 
               GAP=1e+07, yaxis=NULL,
               col2plot="ModelFreq", ...)
}
###

pdf("~/Documents/Heterosis_GWAS/HGWAS_proj/reports/bayes_6traits.pdf", width=10, height=5)
qplot2(pval=cw2, main="CW")
qplot2(pval=cd2, main="CD")
qplot2(pval=cl2, main="CL")
qplot2(pval=akw2, main="AKW")
qplot2(pval=tkw2, main="TKW")
qplot2(pval=kc2, main="KC")
dev.off()


