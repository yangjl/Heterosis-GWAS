### Jinliang Yang
### Joost PNAS Data



tav_map2_jeff <- function(postav, map){
    
    tav <- postav[, 1:2]
    tav$chr <- gsub("_.*", "", tav$snpid)
    tav$chr <- as.numeric(as.character(gsub("S", "", tav$chr)))
    tav$pos <- as.numeric(as.character(gsub("S.*_", "", tav$snpid)))
    
    tav$snp50id <- "A"
    tav$snp50pos <- -9
    
    map$SNP <- as.character(map$SNP)
    
    for(i in 1:nrow(tav)){
        sub <- subset(map, chr == tav$chr[i])
        sub$diff <- abs(sub$pos - tav$pos[i])
        idx <- which.min(sub$diff)
        tav$snp50id[i] <- map[idx, ]$SNP
        tav$snp50pos[i] <- map[idx, ]$pos
        
    }
    tav$diff <- abs(tav$snp50pos - tav$pos)
    return(tav)
}
#################################
jst <- read.delim("largedata/Joostdata/genotype_data")
map <- read.delim("largedata/Joostdata/SNP_data")
grp <- read.csv("largedata/Joostdata/sd01_groups.csv")


tav <- read.csv("cache/TAV_filtered_4plot.csv")
#tav <- read.csv("largedata/TAV_filtered.csv")
postav <- subset(tav, Effect_D > 0)

newtav <- tav_map2_jeff(postav, map)
nntav <- subset(newtav, diff > 1000000)
nntav$snp50id <- gsub("-", ".", nntav$snp50id)


recode_geno <- function( g=subset(grp, GROUP==1 | GROUP==2)){
    geno <- jst[, nntav$snp50id]
    gid <- row.names(geno)
    geno <- apply(geno, 2, as.character)
    geno[geno=="AA"] <- 0
    geno[geno=="AB"] <- 1
    geno[geno=="BB"] <- 2
    geno[geno=="--"] <- 3
    row.names(geno) <- gid
    
    e0g1 <- geno[subset(g, era ==0 & GROUP==1)$GENO_ID,]
    e1g1 <- geno[subset(g, era ==1 & GROUP==1)$GENO_ID,]
    e2g1 <- geno[subset(g, era ==2 & GROUP==1)$GENO_ID,]
    e3g1 <- geno[subset(g, era ==3 & GROUP==1)$GENO_ID,]
    
    e0g2 <- geno[subset(g, era ==0 & GROUP==2)$GENO_ID,]
    e1g2 <- geno[subset(g, era ==1 & GROUP==2)$GENO_ID,]
    e2g2 <- geno[subset(g, era ==2 & GROUP==2)$GENO_ID,]
    e3g2 <- geno[subset(g, era ==3 & GROUP==2)$GENO_ID,]
    
    #frq10 <- getfrq(df = e0g1) #0.5556962
    frq11 <- getfrq(df = e1g1) #0.5516906
    frq12 <- getfrq(df = e2g1) # 0.5579882
    frq13 <- getfrq(df = e3g1) #0.5592332
    
    frq20 <- getfrq(df = e0g2) #0.5595461
    frq21 <- getfrq(df = e1g2) #0.5639241
    frq22 <- getfrq(df = e2g2) #0.5532783
    frq23 <- getfrq(df = e3g2) #0.5557792
    
    
    boxplot(frq10$f)
}

getfrq <- function(df = subgeno1){
    mx <- apply(df, 2, as.numeric)
    f <- apply(mx, 2, function(x) {
        x <- x[x!=3]
        return(sum(x)/(2*length(x)))
    })
    return(as.data.frame(f))
}
