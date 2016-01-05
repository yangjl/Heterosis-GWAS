# purpose: merge and vis the Diallel pheno
# Jinliang Yang
# 2.1.2012

setwd("/Users/yangjl/Documents/Heterosis_GWAS/pheno2011");
f1 <- read.csv("pheno_diallel_master_BLUE.csv");
parents <- read.csv("pheno_namparents_BLUE.csv")
heterosis <- read.csv("pheno_diallel_master_BLUE_heterosis.csv")
sca <- read.csv("pheno_diallel_master_BLUP_SCA.csv")
gca <- read.csv("pheno_diallel_master_BLUP_GCA.csv")
f1$type <- "F1"
parents$type <- "P"

pheno <- rbind(f1, parents)
### Normalize the pheno ###
for (i in 2:8){
  big <- max(pheno[, i])
  small <- min(pheno[,i])
  pheno[, i] <- (pheno[,i]-small)/(big-small);
}


####################Trait per se ##########################################
par(mfrow=c(2,4))
cols <- c("red", "bisque4","blue3","cadetblue4","chartreuse4","chocolate4", "darkslategray")
for (i in 2:8){
  plot(density(subset(pheno, type=="F1")[,i]), col=cols[i-1], xlab="Normalized Scale", lwd=3, xlim=c(-0.3, 1.3), ylim=c(0, 8), main=names(pheno)[i])
  lines(density(subset(pheno, type=="P")[,i]), col=cols[i-1], lty=2, lwd=3)
  #abline(v=subset(pheno, Genotype=="B73")[,i], col=cols[i-2], lwd=3) 
}

### Plot the normalized density plot ###
plot(density(subset(pheno, type=="F1")$KRN), main="All traits", xlab="Normalized Scale", col="red", lwd=2, xlim=c(-0.3, 1.3), ylim=c(0, 8))
lines(density(subset(pheno, type=="P")$KRN), col="red", lty=2, lwd=2)

for (i in 3:8){
  lines(density(subset(pheno, type=="F1")[,i]), col=cols[i-1], lwd=2)
  lines(density(subset(pheno, type=="P")[,i]), col=cols[i-1], lty=2, lwd=2)
  #abline(v=subset(pheno, Genotype=="B73")[,i], col=cols[i-2], lwd=3) 
}
#abline(v=subset(pheno, Genotype=="B73")$KRN, col="red", lwd=3)
###################### HPH & MPH ########################################
### Normalizing the heterosis ###
for (i in 6:26){
  big <- max(heterosis[, i])
  small <- min(heterosis[,i])
  heterosis[, i] <- (heterosis[,i]-0)/(big-small);
}

par(mfrow=c(2,4))
cols <- c("red", "bisque4","blue3","cadetblue4","chartreuse4","chocolate4", "darkslategray")
hph <- heterosis[, c(13,15,17,19,21,23,25)]
mph <- heterosis[, c(14,16,18,20,22,24,26)]
nm <- c("KRN HPH & MPH", "CL HPH & MPH", "CW HPH & MPH", "CD HPH & MPH", "AKW HPH & MPH", "TKW HPH & MPH", "KC HPH & MPH")
for (i in 1:7){
  plot(density(mph[,i]), col=cols[i], lty=2, lwd=3, main=nm[i])
  lines(density(hph[,i]), col=cols[i], lwd=3)
  #abline(v=subset(pheno, Genotype=="B73")[,i], col=cols[i-2], lwd=3) 
}

source("~/Documents/Codes/reusableRcode.r")
pairs(heterosis[, 6:26], text.panel = diag, upper.panel=panel.smooth, 
      lower.panel=panel.cor, gap=0, main="Correlation Plots of the HPH")

####################### GCA ##############################################################################
### Normalizing the heterosis ###
for (i in 2:8){
  big <- max(gca[, i])
  small <- min(gca[,i])
  gca[, i] <- (gca[,i]-small)/(big-small);
}

gca$mean <- apply(gca[,2:8], 1, mean)
gca <- merge(namid, gca, by.x="pop", by.y="Genotype")
gca <- gca[order(gca$mean, decreasing=TRUE),]

#######
#cols <- c("greenyellow", "honeydew3","indianred1", "khaki2", "lightblue3","mediumorchid3", "mediumorchid1")
cols <- c("lightcoral", "bisque4","royalblue3","plum","chartreuse4","chocolate4", "yellow3")
par(mfrow=c(1,1), las=1,col.axis="royalblue")
hh <- t(gca[,3:9]);
mybarcol <- "gray20"
mp <- barplot(hh, beside = TRUE,
              col = cols,
              names.arg = gca$parent,
              legend = rownames(hh), ylim= c(0,1),
              main = "GCA of NAM parents", font.main = 4,
              col.sub = mybarcol,
              cex.names = 0.8)

#segments(mp, hh, mp, hh + 2*sqrt(1000*hh/100), col = mybarcol, lwd = 1.5)
#stopifnot(dim(mp) == dim(hh))# corresponding matrices
mtext(side = 1, at = colMeans(mp), line = 2,
      text = round(colMeans(hh),2), col = "red")

####################### SCA ##############################################################################
### Normalizing the heterosis ###
for (i in 2:8){
  big <- max(sca[, i])
  small <- min(sca[,i])
  sca[, i] <- (sca[,i]-0)/(big-small);
}

#######
par(mfrow=c(2,4))
#cols <- c("greenyellow", "honeydew3","indianred1", "khaki2", "lightblue3","mediumorchid3", "mediumorchid1")
cols <- c("lightcoral", "bisque4","royalblue3","plum","chartreuse4","chocolate4", "yellow3")
for (i in 2:8){
  plot(density(sca[,i]), col=cols[i-1], lty=2, lwd=3, main=names(sca)[i], xlim=c(-1, 1),ylim=c(0,4), xlab="normalized scale")
  lines(density(sca[,i]), col=cols[i-1], lwd=3)
  #abline(v=subset(pheno, Genotype=="B73")[,i], col=cols[i-2], lwd=3) 
}

################################################################

