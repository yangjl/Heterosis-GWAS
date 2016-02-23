### Jinliang Yang
## 


files <- list.files(path="/home/NSF-SAM-GS/HGWAS/SNPdiallel/drealrun", pattern="out1", full.names=TRUE)


source("~/Documents/Github/zmSNPtools/Rcodes/collect_gsout.R")

res <- collect_gsout(dir = "/home/NSF-SAM-GS/HGWAS/SNPdiallel/drealrun", fileptn ="out")
res$trait <- gsub("_.*", "", res$file)
res$effect <- gsub("_run.*", "", res$file)
res$effect <- gsub(".*_", "", res$effect)
write.table(res, "cache/diallel_h2g.csv", sep=",", row.names=FALSE, quote=FALSE)

######
h2g <- read.csv("cache/diallel_h2g.csv")
h2g$trait <- toupper(h2g$trait)
htable <- read.csv("cache/htable_diallel.csv")


h2g_add <- subset(h2g, effect == "add")
h2g_add <- merge(h2g_add, htable, by="trait")


h2g_dom <- subset(h2g, effect == "dom")
h2g_dom <- merge(h2g_dom, htable, by="trait")



lm1 <- lm(h2g_add$h2~h2g_add$pHPH)
p_conf1 <- predict(lm1,interval="confidence")
#p_pred1 <- predict(lm1,interval="prediction")
plot(h2g_add$pHPH, h2g_add$h2, ylim=c(0, 1))
abline(lm1, lwd=2, col="red") ## fit
matlines(sort(h2g_add$pHPH),p_conf1[,c("lwr","upr")], col="red", lty=2, lwd=2, type="b", pch="+")
lines(lowess(h2g_add$pHPH, h2g_add$h2))


lw1 = loess(h2 ~ pHPH, data=h2g_add, degree = 10)
plot(h2 ~ pHPH, data=h2g_add, pch=19, cex=0.1, ylim=c(0, 1))
text(h2g_add$pHPH, h2g_add$h2, labels= h2g_add$trait, cex= 0.7, pos=3)
lines(sort(h2g_add$pHPH), lw1$fitted,col="blue",lwd=3)


points(h2g_dom$pHPH, h2g_dom$h2 - h2g_add$h2, col="red")
points(h2g_dom$pHPH, h2g_dom$h2, col="green")


cor.test(h2g_add$pHPH, h2g_add$h2)
cor.test(h2g_dom$pMPH, h2g_dom$h2 - )

cor.test(h2g_dom[-3, ]$pHPH, h2g_dom[-3, ]$h2 - h2g_add[-3, ]$h2)
cor.test(h2g_dom$pHPH, h2g_dom$h2 - h2g_add$h2)

