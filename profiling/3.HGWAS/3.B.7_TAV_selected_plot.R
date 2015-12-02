### Jinliang Yang
### select the trait-associated variants using a thinning procedure
### 8/15/2014

######
tav <- read.csv("~/Documents/Heterosis_GWAS/HGWAS_proj/reports/StableN.tav758.csv")


pdf("~/Documents/Heterosis_GWAS/HGWAS_proj/graphs/F2_catalog.pdf", height=10, width=10)
# start data plotting:
# plot the frame of the figure
cl <- read.csv("~/Documents/Rcodes/chr_length_B73v2.csv")
plot(c(0, max(cl$BP)/1000000), c(0,105), type="n", main="", 
     xlab="Physical position (Mb)", ylab="", yaxt="n", bty="n")

#### chromosome
for (i in 1:10){
  lines(c(0, cl[i,]$BP/1000000),c(105-10*(i-1), 105-10*(i-1)), lwd=4, col="grey")
  #lines (c(centromere[i,]$Start,centromere[i,]$End),
  #       c(105-10*i, 105-10*i),lwd=5, col="tomato") 
}

axis(side=2, tick =FALSE, las=1, at=c(105, 95, 85, 75, 65, 55, 45, 35, 25, 15), 
     labels=paste("Chr", 1:10, sep=""))

trait <- c("KRN", "CD", "CL", "CW", "AKW", "TKW", "KC")
cols <- c("deeppink4", "burlywood4", "cadetblue4", "chocolate4",
          "cornflowerblue", "darkorchid2", "yellowgreen")
for (t in 1:7){
  myt <- trait[t]
  sub1 <- subset(tav, trait == myt)
  sub1 <- sub1[!duplicated(sub1$bin),]
  for(i in 1:10){
    mysub <- subset(sub1, chr==i)
    if(nrow(mysub) >0){
      points(mysub$pos/1000000, rep(104 -t -10*(i-1), nrow(mysub)), 
             pch=19, cex=1, col=cols[t])
    }
    
  }
}

legend("bottomright", pch=19, cex=1.2, col=cols, legend=trait)
dev.off()

