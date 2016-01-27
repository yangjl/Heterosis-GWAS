### Jinliang Yang
### Nov 25th, 2014
### purpose: cv with recombination rate suggested by reviewer two


ob <- load("cache/cvsnp_p2g_TAV758.RData")

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
cvsnp <- cvsnp[!duplicated(cvsnp$marker),]


####### chisq
M <- as.table(rbind(c(13, 33), c(33, 24)))
dimnames(M) <- list(type = c("big1","small1"),
                    val = c("CV","notCV"))
p1 <- chisq.test(M, corr=FALSE)

#http://bib.oxfordjournals.org/content/9/1/1.full.pdf+html
#gc <- c(10, 190, 800, 3, 100, 900)
#ac <- c(2*gc[1] +gc[2], gc[2]+2*gc[3], 2*gc[4]+gc[5], gc[5]+2*gc[6])
#p1 <- chisq.test(matrix(ac, ncol=2, byrow=T), corr=FALSE)$statistic

eff <- read.csv("cache/TAV_add_dom_eff.csv")
cvsnp <- merge(cvsnp, eff, by.x="marker", by.y="snpid")
######################
pdf("graphs/Test_dominant_recomb.pdf", width=16, height=8)
par(mfrow=c(2,5))
for(i in 1:10){
    submap <- subset(nmap, Chr == i)
    plot(x=submap$Physical/1000000, y=submap$Genetic, pch=16, col="grey", xlab="physical", 
         ylab="genetic", main=paste("Chr", i, sep=""), bty="n")
    
    #### all dots
    subsnp <- subset(cvsnp, chr == i & k > 0.5)
    if(nrow(subsnp) > 0){
        points(x=subsnp$physical/1000000, y=subsnp$genetic, pch=16, col="blue", cex=1.2)
    }
}
dev.off()









