### Jinliang Yang
### stepwise
### 8/13/2014

readSW <- function(infile="AKW_v2_stepwise.mrkRes1"){

  pval <- read.table(infile, header=TRUE)
  pval$chr <- as.numeric(as.character(gsub("_.*", "", pval$ID))) 
  pval$pos <- as.numeric(as.character(gsub(".*_", "", pval$ID)))
  
  message(sprintf("TOTAL SNPs [ %s ]", nrow(pval)))
  return(pval)
}

########
setwd("~/Documents/Heterosis_GWAS/Method/Stepwise/")

krn3 <- readSW(infile="KRN_sw_v3.mrkRes1")
tkw3 <- readSW(infile="TKW_v2_stepwise.mrkRes1")
akw3 <- readSW(infile="AKW_v2_stepwise.mrkRes1")
kc3 <- readSW(infile="KC_v2_stepwise.mrkRes1")
cd3 <- readSW(infile="CD_v2_stepwise.mrkRes3")
cl3 <- readSW(infile="CL_v2_stepwise.mrkRes1")
cw3 <- readSW(infile="CW_v2_stepwise.mrkRes1")

source("~/Documents/Rcodes/save.append.R")
save.append(list=c("tkw3", "akw3", "kc3", "cl3", "cd3", "cw3", "krn3"), 
            file="~/Documents/Heterosis_GWAS/HGWAS_proj/cache/swpval.RData",
            description="")


#### quick check the SNPTEST results
source("~/Documents/Rcodes/quickMHTplot.R")

qplot3 <- function(pval=krn2, ...){
  
  message(sprintf("# of ModelFreq > 0.0001 [ %s ]", nrow(pval)))
  quickMHTplot(res=pval, cex=.6, pch=19, col=rep(c("slateblue", "cyan4"), 5), 
               GAP=1e+07, yaxis=NULL,
               col2plot="tvalue", ...)
}
###

ob3 <- load("~/Documents/Heterosis_GWAS/HGWAS_proj/cache/swpval.RData")

pdf("~/Documents/Heterosis_GWAS/HGWAS_proj/reports/sw_7traits.pdf", width=10, height=5)
qplot3(pval=krn3, main="KRN")
qplot3(pval=cw3, main="CW")
qplot3(pval=cd3, main="CD")
qplot3(pval=cl3, main="CL")
qplot3(pval=akw3, main="AKW")
qplot3(pval=tkw3, main="TKW")
qplot3(pval=kc3, main="KC")
dev.off()


