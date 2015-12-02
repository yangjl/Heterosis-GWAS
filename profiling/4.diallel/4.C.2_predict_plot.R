# Jinliang Yang
# updated: 7/14/2014
# Purpose: 5-fold cross-validation

####################
SAM_predict <- function(ghat="/home/NSF-SAM-GS/GenSel/predictSAM/SAM_predict_large.ghat2",
                        cgrRes="/home/NSF-SAM-GS/GenSel/SAM_run41000_pi.cgrRes1", ...){
  h <- read.table(ghat, skip=1)
  #PEV=var(g/y)
  names(h) <- c("Genotype", "gHat", "log2vol", "Fix")
  p <- read.table(cgrRes, skip=2)
  h$gHat <- h$gHat + p$V2
   
  hist(h$gHat, breaks=30, ...)
  return(h[order(h$gHat),])
}

### 
samlarge <- SAM_predict(ghat="/home/NSF-SAM-GS/GenSel/predictSAM/SAM_predict_large.ghat2",
                        cgrRes="/home/NSF-SAM-GS/GenSel/SAM_run41000_pi.cgrRes1", 
                        xlim=c(19.4,20.6), main="SAM large")

samsmall <- SAM_predict(ghat="/home/NSF-SAM-GS/GenSel/predictSAM/SAM_predict_small.ghat2",
                        cgrRes="/home/NSF-SAM-GS/GenSel/SAM_run41000_pi.cgrRes1", 
                        xlim=c(19.4,20.6), main="SAM small")

################
par(mfrow=c(1,2))
hist(samlarge$gHat, breaks=40, main="SAM large", xlab="log2(volumn)",
     xlim=c(19.4,20.6), col="lightgrey")
abline(v=mean(samlarge$gHat), lwd=2, col="red", lty=2)

hist(samsmall$gHat, breaks=30, main="SAM small", xlab="log2(volumn)",
     xlim=c(19.4,20.6), col="lightgrey")
abline(v=mean(samsmall$gHat), lwd=2, col="red", lty=2)

samlarge$P1 <- gsub("x.*", "", samlarge$Genotype)
samlarge$P2 <- gsub(".*x", "", samlarge$Genotype)
write.table(samlarge, "~/Documents/SAM_GS/SAM_proj/reports/SAM_GS_large.csv", sep=",",
            row.names=FALSE, quote=FALSE)

samsmall$P1 <- gsub("x.*", "", samsmall$Genotype)
samsmall$P2 <- gsub(".*x", "", samsmall$Genotype)
write.table(samsmall, "~/Documents/SAM_GS/SAM_proj/reports/SAM_GS_small.csv", sep=",",
            row.names=FALSE, quote=FALSE)
