### Jinliang Yang
### select the trait-associated variants using a thinning procedure
### 8/15/2014

######
ob <- load(file="~/Documents/Heterosis_GWAS/HGWAS_proj/cache/TAVs_7traits_v2.RData")

## format: snpid, chr, pos, value, bin,  method and trait
nms <- c("snpid", "chr", "pos", "value", "bin", "method")
krn1 <- lskrn[[1]][, c("snpid", "chr", "pos", "log10p", "bin", "Method")]
krn2 <- lskrn[[2]][, c("snpid", "chr", "pos", "ModelFreq", "bin", "Method")]
krn3 <- lskrn[[3]][, c("snpid", "chr", "pos", "Pvalue", "bin", "Method")]
names(krn1) <- nms
names(krn2) <- nms
names(krn3) <- nms
krn0 <- rbind(krn1, krn2, krn3)
krn0$trait <- "KRN"

#####################
FormatOther <- function(myls=lscd, trait="CD"){
  nms <- c("snpid", "chr", "pos", "value", "bin", "method")
  t1 <- myls[[1]][, c("snpid", "chr", "pos", "log10p", "bin")]
  t1$method <- "single-variant"
  names(t1) <- nms
  t2 <- myls[[2]][, c("ID", "chr", "pos", "ModelFreq", "bin")]
  t2$method <- "bayes"
  names(t2) <- nms
  t3 <- myls[[3]][, c("ID", "chr", "pos", "Pvalue", "bin")]
  t3$method <- "stepwise"
  names(t3) <- nms
  
  myt <- rbind(t1, t2, t3)
  myt$trait <- trait
  message(sprintf("snptest [%s], bayes [%s] and stepwise [%s]", 
                  nrow(t1), nrow(t2), nrow(t3)))
  return(myt)
} 

######
cd0 <- FormatOther(myls=lscd, trait="CD")
cl0 <- FormatOther(myls=lscl, trait="CL")
cw0 <- FormatOther(myls=lscw, trait="CW")
akw0 <- FormatOther(myls=lsakw, trait="AKW")
tkw0 <- FormatOther(myls=lstkw, trait="TKW")
kc0 <- FormatOther(myls=lskc, trait="KC")


tav <- rbind(krn0, cd0, cl0, cw0, akw0, tkw0, kc0)
length(unique(tav$snpid)) #[1] 758

tav$bin <- paste(tav$chr, round(tav$pos/1000000,0), sep="_")
length(unique(tav$bin)) #524

write.table(tav, "~/Documents/Heterosis_GWAS/HGWAS_proj/reports/StableN.tav758.csv",
            sep=",", row.names=FALSE, quote=FALSE)




