# Jinliang Yang
# July 18th, 2014
# purpose: plot of heterosis

ob <- load("cache/heterosis_traits.RData")

### changed to dialle only
dkrn <- subset(krn, pop == "Diallel")
dcd <- subset(cd, pop == "Diallel")
dcl <- subset(cl, pop == "Diallel")
dcw <- subset(cw, pop == "Diallel")
dakw <- subset(akw, pop == "Diallel")
dkc <- subset(kc, pop == "Diallel")
dtkw <- subset(tkw, pop == "Diallel")


tot <- rbind(krn, cd, cl, cw, akw, kc, tkw)
sum(tot$pHPH >0)/nrow(tot)


#################
#ph <- c(mean(krn$pHPH), mean(cd$pHPH), mean(akw$pHPH), mean(cl$pHPH), mean(cw$pHPH),
#        mean(kc$pHPH), mean(tkw$pHPH))
#pm <- c(mean(krn$pMPH), mean(cd$pMPH), mean(akw$pMPH), mean(cl$pMPH), mean(cw$pMPH),
#        mean(kc$pMPH), mean(tkw$pMPH))

#################
ph <- c(median(dkrn$pHPH), median(dcd$pHPH), median(dakw$pHPH), median(dcl$pHPH), median(dcw$pHPH),
        median(dkc$pHPH), median(dtkw$pHPH))
ph_mean <- c(mean(dkrn$pHPH), mean(dcd$pHPH), mean(dakw$pHPH), mean(dcl$pHPH), mean(dcw$pHPH),
        mean(dkc$pHPH), mean(dtkw$pHPH))
pm <- c(median(dkrn$pMPH), median(dcd$pMPH), median(dakw$pMPH), median(dcl$pMPH), median(dcw$pMPH),
        median(dkc$pMPH), median(dtkw$pMPH))

htable <- data.frame(trait=c("KRN", "CD", "AKW", "CL", "CW", "KC", "TKW"),
                     pHPH=ph, pHPH_mean=ph_mean, pMPH=pm)


htable <- htable[order(htable$pHPH),]
### mean of the everything
htable$pHPH <- round(htable$pHPH*100, 1)
htable$pMPH <- round(htable$pMPH*100, 1)
write.table(htable, "cache/htable_diallel.csv", sep=",", row.names=FALSE, quote=FALSE)

htable2 <- rbind(krn, cd, cl, cw, akw, kc, tkw)


htable3 <- rbind(dkrn, dcd, dcl, dcw, dakw, dkc, dtkw)
htable3 <- subset(htable3, pHPH<10)
write.table(htable3, "cache/pHPH_diallel_7traits.csv", sep=",", row.names=FALSE, quote=FALSE)

##################################################################################
pdf("graphs/F1_doh.pdf", width=5, height=5)
plot(1:7, htable3$pMPH, bty="n", xaxt="n", xlab="", ylim=c(0, 500),
     ylab="Degress of Heterosis", pch=22, cex=1, bg="lightgrey")
points(1:7, htable3$pHPH, pch=23, cex=1, bg="lightgrey")
axis(side=1, at=1:7, labels=htable$trait)
legend("topleft", legend =c("MPH", "HPH"), pch=c(22,23), bty="n")
dev.off()

pdf("graphs/F1_doh_v3.pdf", width=5, height=5)
htb1 <- subset(htable3, pHPH<10)
bymed2 <- with(htb1, reorder(trait, pHPH, median))
boxplot(pHPH*100 ~ bymed2, data= htb1, notch=TRUE, 
        col="gold",
        main="", ylab="HPH (100%)", xlab="")
dev.off()



pdf("graphs/SF1_pMPH.pdf", width=6, height=6)
htb2 <- subset(htable3, pHPH<10)
bymed2 <- with(htb2, reorder(trait, pMPH, median))
boxplot(pHPH*100 ~ bymed2, data= htb1, notch=TRUE, 
        col="darkgreen",
        main="", ylab="MPH (100%)", xlab="")
dev.off()

boxplot(len~supp*dose, data=ToothGrowth, notch=TRUE, 
        col=(c("gold","darkgreen")),
        main="Tooth Growth", xlab="Suppliment and Dose")
