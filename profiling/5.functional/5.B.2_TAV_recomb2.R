### Jinliang Yang
### Nov 25th, 2014
### purpose: cv with recombination rate suggested by reviewer two

setwd("~/Documents/KRN_GWAS_v3/GWAS3_proj/")
ob <- load("cache/cvsnp_p2g.RData")

########################
# compute recombination rate per 10cM
RecombRate <- function(nmap=map){
    nmap$bin <- paste(nmap$Chr, round(nmap$Genetic/10, 0) , sep="_")
    
    res <- data.frame()
    for(i in unique(nmap$bin)){
        subnmap <- subset(nmap, bin == i)
        subnmap$rate <- (subnmap$Genetic[nrow(subnmap)] - subnmap$Genetic[1])/(subnmap$Physical[nrow(subnmap)] - subnmap$Physical[1]) *1000000
        res <- rbind(res, subnmap)
    }
    return(res)
}

nmap <- RecombRate(nmap=map) #cM/Mb in every 10cM interval
hist(nmap$rate, breaks=50)
nrow(subset(nmap, rate < 1)) #633 
nrow(subset(nmap, rate > 1)) #376

cmrate <- nmap[!duplicated(nmap$bin), ]



#######################
cvsnp$bin <- paste(cvsnp$chr, round(cvsnp$genetic/10, 0) , sep="_")
cvsnp <- merge(cvsnp, cmrate, by="bin", all.x=TRUE)

nrow(subset(cvsnp, type != "Control"))

cvsnp2 <- subset(cvsnp, !is.na(rate) & type != "Control") #103
cvyes <- subset(cvsnp2, cvd==1)
nrow(subset(cvyes, rate >= 1)) #13
nrow(subset(cvyes, rate < 1)) #33

cvno <- subset(cvsnp2, cvd!=1)
nrow(subset(cvno, rate >= 1)) #33
nrow(subset(cvno, rate < 1)) #24

####### chisq
M <- as.table(rbind(c(13, 33), c(33, 24)))
dimnames(M) <- list(type = c("big1","small1"),
                    val = c("CV","notCV"))
p1 <- chisq.test(M, corr=FALSE)

#http://bib.oxfordjournals.org/content/9/1/1.full.pdf+html
#gc <- c(10, 190, 800, 3, 100, 900)
#ac <- c(2*gc[1] +gc[2], gc[2]+2*gc[3], 2*gc[4]+gc[5], gc[5]+2*gc[6])
#p1 <- chisq.test(matrix(ac, ncol=2, byrow=T), corr=FALSE)$statistic

######################
pdf("graphs/cv_recomb.pdf", width=16, height=8)
par(mfrow=c(2,5))
for(i in 1:10){
    submap <- subset(nmap, Chr == i)
    plot(x=submap$Physical/1000000, y=submap$Genetic, pch=16, col="grey", xlab="physical", 
         ylab="genetic", main=paste("Chr", i, sep=""), bty="n")
    
    #### all dots
    subcv <- subset(cvsnp, chr == i)
    if(nrow(subcv) > 0){
        points(x=subcv$physical/1000000, y=subcv$genetic, pch=16, col="blue", cex=1.2)
    }
    #### cvd dots
    subcvd <- subset(subcv, cvd==1)
    if(nrow(subcvd) > 0){
        points(x=subcvd$physical/1000000, y=subcvd$genetic, pch=16, col="red", cex=1.2)
    }
}
dev.off()









