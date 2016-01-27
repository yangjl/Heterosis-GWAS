### Jinliang Yang
### 8/15/2014

#ob <- load("largedata/lcache/TAVs_7traits.RData")
ob <- load("largedata/lcache/TAVs_7traits_v2.RData")
getsnpid <- function(ob=lstkw){
  
  names(ob[[3]])[1] <- "snpid"
  names(ob[[2]])[1] <- "snpid"
  snpid <- c(ob[[1]]$snpid, ob[[2]]$snpid, ob[[3]]$snpid)
  message(sprintf("Total: [%s]; Unique: [%s]", length(snpid), length(unique(snpid))))
  return(snpid)
}
############################
snpid1 <- getsnpid(ob=lstkw)
snpid2 <- getsnpid(ob=lsakw)
snpid3 <- getsnpid(ob=lskc)
snpid4 <- getsnpid(ob=lscd)
snpid5 <- getsnpid(ob=lscw)
snpid6 <- getsnpid(ob=lscl)
snpid7 <- getsnpid(ob=lskrn)

allsnpid <- c(snpid1, snpid2, snpid3, snpid4, snpid5, snpid6, snpid7);
allsnpid <- unique(allsnpid)
#758
write.table(t(allsnpid), "largedata/SNP/TAV_snpid.txt",
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")








