## thinning the GWAS results

TAVselect <- function(snptest=cw, bayes=cw2, sw=cw3, 
                      bayescutoff=0.02, snptestcutoff=20, cwcutoff=2,
                      binsize=1e7, topinbin=10){
  
  message(sprintf("# input: SNPTEST [%s], Bayes [%s] and SW [%s]",
                  nrow(snptest), nrow(bayes), nrow(sw) ))
  #####CUTOFFs
  try(snptest$log10p <- -log10(snptest$pval))
  snptest <- subset(snptest, log10p > snptestcutoff)
  
  bayes <- subset(bayes, ModelFreq > bayescutoff)
  
  sw$log10p <- -log10(sw$Pvalue)
  sw <- subset(sw, log10p > cwcutoff)
  
  message(sprintf("# >cutoff: SNPTEST [%s], Bayes [%s] and SW [%s]",
                  nrow(snptest), nrow(bayes), nrow(sw) ))
  #sw <- subset(sw, )
  snptest$bin <- paste(snptest$chr, round(snptest$pos/binsize, 0), sep="_")
  bayes$pos <- as.numeric(as.character(bayes$pos))
  bayes$bin <- paste(bayes$chr, round(bayes$pos/binsize, 0), sep="_")
  sw$bin <- paste(sw$chr, round(sw$pos/binsize, 0), sep="_")
  
  ### THINNING
  snptest <- snptest[order(snptest$log10p, decreasing=TRUE),]
  out1 <- data.frame()
  for(bini in unique(snptest$bin)){
    tem1 <- subset(snptest, bin==bini)
    if(nrow(tem1) > topinbin){
      tem1 <- tem1[1:topinbin, ]
    }
    out1 <- rbind(out1, tem1)
  }
  ### Bayes
  bayes <- bayes[order(bayes$ModelFreq, decreasing=TRUE),]
  out2 <- data.frame()
  for(bini in unique(bayes$bin)){
    tem2 <- subset(bayes, bin==bini)
    if(nrow(tem2) > topinbin){
      tem2 <- tem2[1:topinbin, ]
    }
    out2 <- rbind(out2, tem2)
  }
  ### stepwise
  #snptest <- snptest[order(snptest$log10p, decreasing=TRUE),]
  if(sum(names(sw) %in% "tvalue") >0){
    sw <- sw[order(sw$tvalue, decreasing=TRUE),]
  }else{
    sw <- sw[order(sw$log10p, decreasing=TRUE),]
  }
  
  out3 <- data.frame()
  for(bini in unique(sw$bin)){
    tem3 <- subset(sw, bin==bini)
    if(nrow(tem3) > topinbin){
      tem3 <- tem3[1:topinbin, ]
    }
    out3 <- rbind(out3, tem3)
  }
  message(sprintf("# after thinning: SNPTEST [%s], Bayes [%s] and SW [%s]",
                  nrow(out1), nrow(out2), nrow(out3) ))
  myout <- list(snptest=out1, bayes=out2, stepwise=out3)
  return(myout)
  
}