### Jinliang Yang
### Feb 23th, 2016
### plot h2g and heterosis


source("~/Documents/Github/zmSNPtools/Rcodes/collect_gsout.R")

res <- collect_gsout(dir = "/home/NSF-SAM-GS/HGWAS/SNPdiallel/drealrun", fileptn ="out")
res$trait <- gsub("_.*", "", res$file)
res$effect <- gsub("_run.*", "", res$file)
res$effect <- gsub(".*_", "", res$effect)
write.table(res, "cache/diallel_h2g.csv", sep=",", row.names=FALSE, quote=FALSE)

######
htable <- read.csv("cache/htable_diallel.csv")

get_h2 <- function(){
    h2g <- read.csv("cache/diallel_h2g.csv")
    h2g$trait <- toupper(h2g$trait)
    
    h2g_add <- subset(h2g, effect == "add")
    h2g_dom <- subset(h2g, effect == "dom")
    
    h2 <- merge(h2g_add[, c("trait", "h2")], h2g_dom[, c("trait", "h2")], by="trait")
    names(h2) <- c("trait", "h2a", "h2")
    h2$h2d <- h2$h2 - h2$h2a
    h2 <- merge(h2, htable, by="trait")
    
    return(h2)
}

#####
h2 <- get_h2()
h2 <- h2[order(h2$pHPH),]
th2 <- t(h2[, c(2,4)])
names(th2) <- h2$trait
barplot(th2, main="Car Distribution by Gears and VS",
        xlab="Number of Gears", col=c("darkblue","red"),
        legend = rownames(th2))

library(ggplot2)
c <- ggplot(h2, aes(pHPH, h2a))
c + stat_smooth()


###################################
h2 <- get_h2()
h2 <- h2[order(h2$pHPH),]

pdf("graphs/h2g_heterosis.pdf", width=5, height=5)
smoothingSpline = smooth.spline(x=h2$pHPH, y=h2$h2a, spar=0.8)
plot(x=h2$pHPH, y=h2$h2a, pch=16, type="n", ylim=c(0,1), xlim=c(-10, 300), xlab="HPH (100%)", ylab="Variance Explained")
text(x=h2$pHPH, y=h2$h2a, labels=h2$trait, col="red", offset=0.6, cex=0.7)
lines(smoothingSpline, lwd=3, col="red")

smoothingSpline2 = smooth.spline(x=h2$pHPH, y=h2$h2d, spar=0.8)
#plot(x=h2$pHPH, y=h2$h2a, pch=16, type="n", ylim=c(0,1))
text(x=h2$pHPH, y=h2$h2d, labels=h2$trait, col="blue", offset=0.6, cex=0.7)
lines(smoothingSpline2, lwd=3, col="blue")

ss3 = smooth.spline(x=h2$pHPH, y=h2$h2, spar=0.8)
#plot(x=h2$pHPH, y=h2$h2a, pch=16, type="n", ylim=c(0,1))
text(x=h2$pHPH, y=h2$h2, labels=h2$trait, col="black", offset=0.6, cex=0.7)
lines(ss3, lwd=3, col="black")
dev.off()








