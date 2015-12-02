## Jinliang Yang
## 7/21/2014
## plot synergistic effect QTLs

#########
synPlot <- function(sqtl=stql, GAP=10000000, trait1="KRN", trait2="AKW",
                    col1="red", col2="blue", ...){
  
  source("~/Documents/Rcodes/newpos.R")
  #res <- newpos(res, GAP = GAP)
  chrtick <- chrline_tick(GAP = GAP)
  namp <- read.csv("~/DBcenter/VariationDB/NAM_populations.csv")
  namp$parent <- toupper(namp$parent)
  namp$ypos <- 26:1
  
  #### change RIL to MO17
  sqtl[sqtl$pop == "RIL", ]$pop <- "MO17" 
  names(sqtl)[4] <- "pos"
  sqtl <- newpos(sqtl, GAP=GAP)
  sqtl <- merge(sqtl, namp[, c("parent", "ypos")], by.x="pop", by.y="parent", all.x=TRUE)
  
  
  sub1 <- subset(sqtl, trait==trait1)
  sub2 <- subset(sqtl, trait==trait2)
  
  sub11 <- subset(sub1, effect > 0)
  sub12 <- subset(sub1, effect < 0)
  sub21 <- subset(sub2, effect > 0)
  sub22 <- subset(sub2, effect < 0)
  
  
  plot(x=-1000, y=-1000,  type="p", xaxt="n", yaxt="n", xlab="", ylab="",
       xlim=c(0, max(chrtick$chrlines)), ylim=c(0, 26.5), bty="n", ...)
  axis(side=2, at=26:1, las=2, labels=namp$parent, tick=FALSE)
  axis(side=1, at=chrtick$ticks, tick=TRUE, labels=c("chr1", "chr2", "chr3", "chr4", "chr5", 
       "chr6", "chr7", "chr8", "chr9", "chr10"))
  abline(v=chrtick$chrlines, col="grey")
  abline(h=1:26, col="grey")
  points(x=sub11$newpos, y=sub11$ypos-0.3, pch=24, col=col1, bg=col1)
  points(x=sub12$newpos, y=sub12$ypos-0.3, pch=25, col=col1, bg=col1)
  points(x=sub21$newpos, y=sub21$ypos-0.7, pch=24, col=col2, bg=col2)
  points(x=sub22$newpos, y=sub22$ypos-0.7, pch=25, col=col2, bg=col2)
  
  legend("topright", inset=c(0, -0.05), c(trait1, trait2, "positive effect", "negative effect"), 
         pch = c(19, 19, 24, 25), x.intersp=0.5, xpd=TRUE,
         horiz = TRUE, cex=0.8, col=c(col1, col2, "black", "black"))
  
}

###########
sqtl <- read.csv("~/Documents/Heterosis_GWAS/HGWAS_proj/reports/sep_qtl_table.csv")
sqtl <- subset(sqtl, pop != "BxRIL" & pop != "MxRIL")

synPlot(sqtl=sqtl, GAP=10000000, trait1="KRN", trait2="AKW",
        col1="red", col2="blue", main="KRN vs. AKW")
synPlot(sqtl=sqtl, GAP=10000000, trait1="KRN", trait2="CL",
        col1="red", col2="darkgoldenrod4", main="KRN vs. CL")
synPlot(sqtl=sqtl, GAP=10000000, trait1="CD", trait2="CL",
        col1="blue3", col2="darkgoldenrod4", main="CD vs. CL")

##################
pdf("~/Documents/Heterosis_GWAS/HGWAS_proj/reports/S.F3_synqtl.pdf", height=7, width=14)
synPlot(sqtl=sqtl, GAP=10000000, trait1="KRN", trait2="AKW",
        col1="red", col2="blue", main="KRN vs. AKW")
synPlot(sqtl=sqtl, GAP=10000000, trait1="KRN", trait2="CL",
        col1="red", col2="darkgoldenrod4", main="KRN vs. CL")
synPlot(sqtl=sqtl, GAP=10000000, trait1="CD", trait2="CL",
        col1="blue3", col2="darkgoldenrod4", main="CD vs. CL")
dev.off()




