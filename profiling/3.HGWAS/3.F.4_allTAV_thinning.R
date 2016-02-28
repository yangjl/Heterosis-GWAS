### Jinliang Yang
### Feb 27th, 2016


################################################################
thinning_TAV <- function(varcutoff=0){
    alltav <- read.csv("largedata/allTAVs_effect.csv")
    alltav$snpid <- paste0("S", alltav$snpid)
    files <- list.files(path="largedata/snpeff", pattern="snpe$", full.names=TRUE)
    output <- data.frame()
    #traits <- c("KRN", "CD", "AKW", "CL", "CW", "KC", "TKW")
    for(i in 1:length(files)){
        h1 <- read.table(files[i], header=TRUE)
        names(h1) <- c("snpid","chr","pos","Effect_A","Effect_D","Effect_A2","Effect_D2","h2_mrk_A", 
                       "h2_mrk_D","H2_mrk","h2_mrk_A_p","h2_mrk_D_p","H2_mrk_p","log10_h2_mrk_A","log10_h2_mrk_D","log10_H2_mrk")
        
        myt <- gsub(".*\\/", "", files[i])
        myt <- gsub("_.*", "", myt)
        subtav <- subset(alltav, trait %in% myt)
        out <- merge(subtav, h1, by="snpid")
        
        out <- out[order(out$H2_mrk, decreasing=TRUE),]
        out <- subset(out, H2_mrk > varcutoff)
        
        out$k <- out$Effect_D/out$Effect_A
        message(sprintf("###>>> processing file: [ %s ] and get [ %s ] SNPs", files[i], nrow(out)))
        output <- rbind(output, out)
    }
    return(output)
}

#### filtered by the median value of the variance explained
TAV <- thinning_TAV(varcutoff=2e-04)
write.table(TAV, "largedata/TAV_filtered.csv", sep=",", row.names=FALSE, quote=FALSE)


