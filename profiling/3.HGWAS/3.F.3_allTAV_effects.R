### Jinliang Yang
### 02/08/2016

###########
getsnpid <- function(myt=lstkw, trait="TKW"){
    
    ## [[1]] SNPTEST
    ## [[2]] BayesC
    ## [[3]] Stepwise
    res1 <- myt[[1]][, c("snpid", "log10p", "beta", "se")]
    names(res1) <- c("snpid", "value", "effect", "var")
    res1$method <- "snptest"
    res2 <- myt[[2]][, c("ID", "ModelFreq","Effect", "EffectVar")]
    names(res2) <- c("snpid", "value", "effect", "var")
    res2$method <- "BayesC"
    res3 <- myt[[3]][, c("ID", "tvalue",  "Effect","Pvalue")]
    names(res3) <- c("snpid", "value", "effect", "var")
    res3$method <- "Stepwise"
    
    res <- rbind(res1, res2, res3)
    res$trait <- trait
    message(sprintf("###>>> Total: [%s]; Unique: [%s]", length(res$snpid), length(unique(res$snpid))))
    return(res)
}
ob <- load(file="largedata/lcache/TAVs_7traits.RData")
snpid1 <- getsnpid(myt=lstkw, trait="TKW")
###>>> Total: [347]; Unique: [337]
snpid2 <- getsnpid(myt=lsakw, trait="AKW")
###>>> Total: [797]; Unique: [780]
snpid3 <- getsnpid(myt=lskc, trait="KC")
###>>> Total: [786]; Unique: [774]
snpid4 <- getsnpid(myt=lscd, trait="CD")
###>>> Total: [1412]; Unique: [1383]
snpid5 <- getsnpid(myt=lscw, trait="CW")
###>>> Total: [893]; Unique: [872]
snpid6 <- getsnpid(myt=lscl, trait="CL")
###>>> Total: [973]; Unique: [962]


names(lskrn[[1]])[4] <- "value"
names(lskrn[[2]])[4] <- "value"
names(lskrn[[3]])[4] <- "value"

snpid7 <- rbind(lskrn[[1]][, c(1,4,5)], lskrn[[2]][, c(1,4,5)], lskrn[[3]][, c(1,4,5)])

snpid7$effect <- NA
snpid7$var <- NA
snpid7$trait <- "KRN"
snpid7 <- snpid7[, c("snpid", "value", "effect", "var", "Method", "trait")]
names(snpid7)[5] <- "method"

allTAVs <- rbind(snpid1, snpid2, snpid3, snpid4, snpid5, snpid6, snpid7)
write.table(allTAVs, "largedata/allTAVs_effect.csv", sep=",", row.names=FALSE, quote=FALSE)

