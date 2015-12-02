# Jinliang Yang
# July 18th, 2014
# purpose: cor plot of 7 phenotypic traits


source("~/Documents/Github/zmSNPtools/Rcodes/Correlation_plot.R")
# read in the data
pheno <- read.table("cache/pheno_matrix.txt", header=TRUE)

pdf("reports/S.F2_cor.pdf", height=8, width=8)
pairs(pheno[, 4:10], text.panel = diag, upper.panel=panel.smooth, 
      lower.panel=panel.cor, gap=0, main="", pch=19, col="grey", lwd=2)
dev.off()
