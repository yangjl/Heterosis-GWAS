### Jinliang Yang
### 01/27/2016



num <- read.csv("cache/TAV_add_dom_num.csv")
eff <- read.csv("cache/TAV_add_dom_eff.csv")

num$file <- gsub(".*\\/|_.*", "", num$file)

library(reshape2)



out2 <- melt(num[, c(c(1, 6,7))], id.var="file")
#############################################
library(ggplot2)
source("~/Documents/Github/zmSNPtools/Rcodes/multiplot.R")

ob <- load("cache/heterosis_traits.RData")

med <- data.frame(trait=c("KRN", "CD", "CL", "CW", "KC", "AKW", "TKW"), 
                   phph= c(mean(krn$pHPH), abs(mean(cd$pHPH)),mean(cl$pHPH), 
                           mean(cw$pHPH),mean(kc$pHPH),mean(akw$pHPH),mean(tkw$pHPH) ))
med <- med[order(med$phph),]







#########################################
out1 <- num[, 1:6]
out1$tot <- apply(out1[, -1], 1, sum)
out1$k1 <- out1$k1/out1$tot
out1$k2 <- out1$k2/out1$tot
out1$k3 <- out1$k3/out1$tot
out1$k4 <- out1$k4/out1$tot
out1$k5 <- out1$k5/out1$tot
out1 <- melt(out1[, 1:6], id.var="file")

out2 <- num[, c(1, 7:9)]
out2$h2_A <- out2$h2_A/out2$h2_tot
out2$h2_D <- out2$h2_D/out2$h2_tot
out2 <- melt(out2[, 1:3], id.var="file")

theme_set(theme_grey(base_size = 18)) 
p1 <- ggplot(out1, aes(x=factor(file, levels=med$trait), y=value, 
                       fill=factor(variable, labels=c("DD > 1", " 0.5 < DD < 1", "-0.5 < DD < 0.5", "-1< DD < -0.5", "DD < -1")))) + 
    geom_bar(position=position_dodge(), stat="identity") +
    xlab("") +
    ylab("Frequency") +
    ggtitle("") + theme_bw() +
    labs(fill="Degree of Dominance") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))

p2 <- ggplot(out2, aes(x=factor(file, levels=med$trait), y=value, 
                       fill=factor(variable, labels=c("Additive", "Dominance")))) + 
    geom_bar(position=position_dodge(), stat="identity") +
    xlab("") +
    ylab("Proportion of Variance Explained") +
    ggtitle("") + theme_bw() +
    labs(fill="Variance Components") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))

#multiplot(p1, p2, cols=2)

########
pdf("graphs/Figure5_DD_var.pdf", width=13, height=5)
multiplot(p1, p2, cols=2)
dev.off()



