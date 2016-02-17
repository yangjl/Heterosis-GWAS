## Jinliang Yang
## July 17th, 2014
## plot the phenotypic traits
# updated: 12.2.2015


plotKRN <- function(krnall){
    krnall$norv <- scale(krnall$value)
    plot(density(subset(krnall, pop=="NAMRIL")$norv), col="black", lwd=4, ylim=c(0, 0.8), bty="n",
         main="Kernel Row Number", xlab="KRN")
    lines(density(subset(krnall, pop=="Diallel")$norv), col="gold", lwd=4)
    lines(density(subset(krnall, pop=="BxRIL")$norv), col="blue", lwd=4)
    lines(density(subset(krnall, pop=="MxRIL")$norv), col="red", lwd=4)
    #abline(v=17.1, col="blue", lty=2, lwd=2)
    temp <- legend("topright", legend =c(" ", " ", " ", " "), title="Population", 
                   text.width=strwidth("MxRIL"),
                   lty=1, lwd=3, col=c("black", "gold", "blue", "red"), 
                   xjust=1, yjust=1)
    text(temp$rect$left + temp$rect$w, temp$text$y,
         c("RIL", "Diallel", "BxRIL", "MxRIL"), pos=2)  
}

###########
denplot <- function(tcol="KRN", ...){
    par(mar=c(2,2,3,2))
    myp <- subset(trait, trait==tcol)
    myp$norv <- scale(myp$value)
    plot(density(subset(myp, pop=="NAMRIL")$norv), col="black", lwd=2, bty="n",
         main=tcol, xlab="", ...)
    lines(density(subset(myp, pop=="Diallel")$norv), col="gold", lwd=2)
    lines(density(subset(myp, pop=="BxRIL")$norv), col="blue", lwd=2)
    lines(density(subset(myp, pop=="MxRIL")$norv), col="red", lwd=2)
}

#######
#nf <- layout(matrix(c(1,1,2,3,4,1,1,5,6,7), 2, 5, byrow = TRUE))
#layout.show(nf)

trait <- read.table("reports/S1_pheno.txt", header=TRUE)

pdf("graphs/S.F1_pheno.pdf", width=10, height=4)
layout(matrix(c(1,1,2,3,4,1,1,5,6,7), 2, 5, byrow = TRUE))
#par(mar=c(0,0,0,0))
plotKRN(krnall=subset(trait, trait=="KRN"))
denplot(tcol="CD", xlim=c(-4, 4), ylim=c(0, 1))
denplot(tcol="AKW", xlim=c(-4, 5), ylim=c(0, 1))
denplot(tcol="CL", xlim=c(-4, 5), ylim=c(0, 0.8))
denplot(tcol="CW", xlim=c(-4, 8), ylim=c(0, 0.8))
denplot(tcol="KC", xlim=c(-4, 6), ylim=c(0, 0.8))
denplot(tcol="TKW", xlim=c(-4, 6), ylim=c(0, 1.2))
dev.off()
