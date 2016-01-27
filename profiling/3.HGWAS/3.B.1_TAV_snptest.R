### Jinliang Yang
###
### 8/13/2014

readSNPTEST <- function(pwd="~/Documents/Heterosis_GWAS/Method/SNPTEST", 
                        snptestfile=c("snptest_TKW_chr","-covars.out"), bychr=TRUE){
  
  setwd(pwd)
  out <- data.frame()
  if(bychr){
    for(i in 1:10){
      infile <- paste(snptestfile[1], i, snptestfile[2], sep="")
      tab5rows <- read.table(infile, nrow=5, header=TRUE)
      classes <- sapply(tab5rows, class)
      
      pval <- read.table(infile, header=TRUE, colClasses=classes)
      pval <- pval[, c(1, 2, 4, 20:23)]
      
      names(pval) <-  c("chr", "snpid", "pos", "pval", "info", "beta", "se")
      out <- rbind(out, pval)
    } 
  }
  message(sprintf("TOTAL SNPs [ %s ]", nrow(out)))
  return(out)
}

########
tkw <- readSNPTEST(pwd="~/Documents/Heterosis_GWAS/Method/SNPTEST", 
                   snptestfile=c("snptest_TKW_chr","-covars.out"), bychr=TRUE)

akw <- readSNPTEST(pwd="~/Documents/Heterosis_GWAS/Method/SNPTEST", 
                   snptestfile=c("snptest_AKW_chr","-covars.out"), bychr=TRUE)

kc <- readSNPTEST(pwd="~/Documents/Heterosis_GWAS/Method/SNPTEST", 
                   snptestfile=c("snptest_KC_chr","-covars.out"), bychr=TRUE)

cl <- readSNPTEST(pwd="~/Documents/Heterosis_GWAS/Method/SNPTEST", 
                   snptestfile=c("snptest_CL_chr","-covars.out"), bychr=TRUE)

cd <- readSNPTEST(pwd="~/Documents/Heterosis_GWAS/Method/SNPTEST", 
                   snptestfile=c("snptest_CD_chr","-covars.out"), bychr=TRUE)

cw <- readSNPTEST(pwd="~/Documents/Heterosis_GWAS/Method/SNPTEST", 
                   snptestfile=c("snptest_CW_chr","-covars.out"), bychr=TRUE)

source("~/Documents/Rcodes/save.append.R")
save.append(list=c("tkw", "akw", "kc", "cl", "cd", "cw"), 
            file="~/Documents/Heterosis_GWAS/HGWAS_proj/cache/snptestpval.RData",
            description="")


#### quick check the SNPTEST results
source("~/Documents/Rcodes/quickMHTplot.R")

qplot1 <- function(pval=tkw, ...){
  pval$log10p <- -log10(pval$pval)
  pval <- subset(pval, log10p > 10)
  message(sprintf("# of -log10p >5 [ %s ]", nrow(pval)))
  quickMHTplot(res=pval, cex=.3, pch=19, col=rep(c("slateblue", "cyan4"), 5), 
               GAP=5e+06, yaxis=NULL,
               col2plot="log10p", ...)
}
###

ob1 <- load("~/Documents/Heterosis_GWAS/HGWAS_proj/cache/snptestpval.RData")
pdf("~/Documents/Heterosis_GWAS/HGWAS_proj/reports/snptest_6traits.pdf", width=10, height=5)
qplot1(pval=cw, main="CW")
# of -log10p >5 [ 158179 ]
qplot1(pval=cd, main="CD")
# of -log10p >5 [ 462525 ]
qplot1(pval=cl, main="CL")
# of -log10p >5 [ 101406 ]
qplot1(pval=akw, main="AKW")
# of -log10p >5 [ 12232 ]
qplot1(pval=tkw, main="TKW")
# of -log10p >5 [ 5189 ]
qplot1(pval=kc, main="KC")
dev.off()


