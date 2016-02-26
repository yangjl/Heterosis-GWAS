# Jinliang Yang
# July 18th, 2014
# purpose: cor plot of 7 phenotypic traits


source("~/Documents/Github/zmSNPtools/Rcodes/Correlation_plot.R")
# read in the data
pheno <- read.table("cache/pheno_matrix.txt", header=TRUE)

pdf("graphs/SF3_cor_traitperse.pdf", height=8, width=8)
pairs(pheno[, c(5,4,6:10)], text.panel = diag, upper.panel=panel.smooth, 
      lower.panel=panel.cor, gap=0, main="", pch=19, col="grey", lwd=2)
dev.off()

pr <- c(0.18, 0.059, 0.71, 0.84, 0.78, 0.9)
mean((pr-r)/r)

htable <- read.csv("cache/htable_diallel.csv")
cor.test(htable$pHPH_mean[-7], pr)


####################################################
hph <- read.csv("cache/pHPH_diallel_7traits.csv")
traits <- c("KRN", "CD", "AKW", "CL", "CW", "KC", "TKW")

p1 <- subset(hph, trait == traits[1])
p1 <- p1[, c("Genotype", "pHPH")]
names(p1)[2] <- traits[1]

for(i in 2:7){
    tem <- subset(hph, trait == traits[i])
    p2 <- tem[, c("Genotype", "pHPH")]
    names(p2)[2] <- traits[i]
    p1 <- merge(p1, p2, by="Genotype")
}

pdf("graphs/SF2_pHPH_corplot.pdf", height=8, width=8)
pairs(p1[, 2:8], text.panel = diag, upper.panel=panel.smooth, 
      lower.panel=panel.cor, gap=0, main="", pch=19, col="grey", lwd=2)
dev.off()

#########################
htable <- read.csv("cache/htable_diallel.csv")
r <- c(0.091, 0.31, 0.4, 0.55, 0.51, 0.82)
cor.test(htable$pHPH_mean[-7], r)


##########################
g1 <- ggplot(p1, aes(x=KRN, y=TKW)) +
    theme_bw() + xlab("KRN") + ylab("TKW") +
    guides(colour=FALSE, linetype=FALSE) +
    geom_point(col="grey", size=0.8) +
    geom_smooth(method="gam", size=1.3) +
    theme(axis.text.y = element_text(angle = 90, hjust = 1))
g2 <- ggplot(p1, aes(x=CD, y=TKW)) +
    theme_bw() + xlab("CD") + ylab("") +
    guides(colour=FALSE, linetype=FALSE) +
    geom_point(col="grey", size=0.8) +
    geom_smooth(method="gam", size=1.3) +
    theme(axis.text.y = element_text(angle = 90, hjust = 1))
g3 <- ggplot(p1, aes(x=AKW, y=TKW)) +
    theme_bw() + xlab("AKW") + ylab("") +
    guides(colour=FALSE, linetype=FALSE) +
    geom_point(col="grey", size=0.8) +
    geom_smooth(method="gam", size=1.3) +
    theme(axis.text.y = element_text(angle = 90, hjust = 1))
g4 <- ggplot(p1, aes(x=CL, y=TKW)) +
    theme_bw() + xlab("CL") + ylab("TKW") +
    guides(colour=FALSE, linetype=FALSE) +
    geom_point(col="grey", size=0.8) +
    geom_smooth(method="gam", size=1.3) +
    theme(axis.text.y = element_text(angle = 90, hjust = 1))
g5 <- ggplot(p1, aes(x=CW, y=TKW)) +
    theme_bw() + xlab("CW") + ylab("") +
    guides(colour=FALSE, linetype=FALSE) +
    geom_point(col="grey", size=0.8) +
    geom_smooth(method="gam", size=1.3) +
    theme(axis.text.y = element_text(angle = 90, hjust = 1))
g6 <- ggplot(p1, aes(x=KC, y=TKW)) +
    theme_bw() + xlab("KC") + ylab("") +
    guides(colour=FALSE, linetype=FALSE) +
    geom_point(col="grey", size=0.8) +
    geom_smooth(method="gam", size=1.3) +
    theme(axis.text.y = element_text(angle = 90, hjust = 1))



#multiplot(p1, p4, p2, p5, p3, p6, cols=3)
source("~/Documents/Github/zmSNPtools/Rcodes/multiplot.R")
pdf("graphs/tem1.pdf", width=10, height=6)
multiplot(g1, g2, g3,g4, g5, g6, cols=3)
dev.off()



