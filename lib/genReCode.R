# Jinliang Yang
# 8/15/2014
# update: sapply

###############
genReCode <- function(gen=sub, code="GenSel", outfile=""){
  
  if(code=="GenSel"){
    aa <- 10;
    ab <- 0;
    an <- 5;
    bb <- -10;
    bn <- -5;
    nn <- 0;
  }else if(code=="PLINK"){
    aa <- "A A";
    ab <- "A B";
    an <- "N N";
    bb <- "B B";
    bn <- "N N";
    nn <- "N N"; 
  }else if(code=="raw"){
    aa <- "AA";
    ab <- "AB";
    an <- "AN";
    bb <- "BB";
    bn <- "BN";
    nn <- "NN"; 
  }
  
  ### create progress bar
  pb <- txtProgressBar(min = 0, max = nrow(gen), style = 3)
  out <- data.frame()
  for(i in 1:nrow(gen)){
    v <- gen$V3[i]
    for(j in 2:(ncol(gen) %/% 3 - 1) ){
      AA = gen[i, j*3+1]
      AB = gen[i, j*3+2]
      BB = gen[i, j*3+3]
      if(AA==1 & AB==0 & BB==0){
        snp <- aa
      }else if(AA==0 & AB==1 & BB==0){
        snp <- ab
      }else if(AA==0.5 & AB==0.5 & BB==0){
        snp <- an
      }else if(AA==0 & AB==0 & BB==1){
        snp <- bb
      }else if(AA==0 & AB==0.5 & BB==0.5){
        snp <- bn
      }else if(AA==0 & AB==0 & BB==0){
        snp <- nn
      }else{
        stop("unexpected pattern!!!")
      }
      v <- c(v, snp)
    }
    vt <- t(v)
    out <- rbind(out, vt)
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  write.table(out, outfile, sep=",", row.names=FALSE, quote=FALSE, col.names=FALSE)
}




