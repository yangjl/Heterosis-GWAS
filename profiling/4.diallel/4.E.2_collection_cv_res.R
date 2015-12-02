### Jinliang Yang
### 9/21/2014
### collect the cv results


setwd("/home/NSF-SAM-GS/HGWAS/SNPdiallel/dcv2/krn")


################
collect_h2 <- function(pwd="/home/NSF-SAM-GS/HGWAS/SNPdiallel/predict/s1_predict",
                       traits=traits, outraw=TRUE){
  ### output file
  out <- data.frame()
  
  for(i in 1:7){
    #mypwd <- paste(pwd, traits[i], sep="")
    files <- list.files(path = pwd, pattern="out1$")    
    
    ###############
    setwd(pwd)
    for(j in 1:length(files)){
      res <- readLines(files[j])
      idx1 <- grep("unweighted correlation", res)
      h2 <- gsub(".* = ", "", res[idx1])
      tem <- data.frame(file=files[j], h2=h2)
      out <- rbind(out, tem)
    }
    
  }
  ##########
  out$trait <- gsub("_.*", "", out$file)
  out$effect <- gsub(".*_", "", out$file)
  ####
  if(outraw){
    return(out) 
  }else{
    out$h2 <- as.numeric(as.character(out$h2))
    myout <- data.frame()
    for(i in 1:length(traits)){
      out2 <- subset(out, trait == traits[i] & effect == "add.out1")
      out3 <- subset(out, trait == traits[i] & effect == "dom.out1")
      tem <- data.frame(trait=traits[i], addh2=mean(out2$h2), addsd=sd(out2$h2),
                        domh2=mean(out3$h2), domsd=sd(out3$h2))
      myout <- rbind(myout, tem)
    }
    return(myout)
  }
}


#############
traits =c("krn", "cd", "akw", "cl", "cw", "kc", "tkw")
cv1 <- collect_h2(pwd="/home/NSF-SAM-GS/HGWAS/SNPdiallel/predict/s1_predict",
                  traits=traits, outraw=FALSE)
cv1$imp <- cv1$domh2 - cv1$addh2

cv2 <- collect_h2(pwd="/home/NSF-SAM-GS/HGWAS/SNPdiallel/predict/s2_predict",
                  traits=traits, outraw=FALSE)
cv2$imp <- cv2$domh2 - cv2$addh2

save(list=c("cv1", "cv2"), file="~/Documents/Heterosis_GWAS/HGWAS_proj/cache/cv_res.RData")
load("~/Documents/Heterosis_GWAS/HGWAS_proj/cache/cv_res.RData")
###########################
par(mfrow=c(1,2))
plot(x=1:7, y=cv1$addh2, type="p", ylim=c(0,1), xlab="", xaxt="n", col="red", lwd=3,
     main="Strategy 1", ylab="Prediction Accuracy (r)")
axis(side=1, at=1:7, labels=toupper(traits) )
lines(x=1:7, y=cv1$domh2, type="p", col="blue", lwd=3)
lines(x=1:7, y=cv1$imp, type="p", col="black", lwd=3, lty=2)

plot(x=1:7, y=cv2$addh2, type="p", ylim=c(0,1), xlab="", xaxt="n", col="red", lwd=3,
     main="Strategy 2", ylab="Prediction Accuracy (r)")
axis(side=1, at=1:7, labels=toupper(traits) )
lines(x=1:7, y=cv2$domh2, type="p", col="blue", lwd=3)
lines(x=1:7, y=cv2$imp, type="p", col="black", lwd=3, lty=2)







