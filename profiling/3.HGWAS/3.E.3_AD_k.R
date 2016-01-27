### Jinliang Yang
### 09/08/2015

get_AD_eff <- function(files, tav){
    
    output1 <- output2 <- data.frame()
    for(i in 1:length(files)){
        h1 <- read.table(files[i], header=TRUE)
        names(h1) <- c("snpid","chr","pos","Effect_A","Effect_D","Effect_A2","Effect_D2","h2_mrk_A", 
                       "h2_mrk_D","H2_mrk","h2_mrk_A_p","h2_mrk_D_p","H2_mrk_p","log10_h2_mrk_A","log10_h2_mrk_D","log10_H2_mrk")
        
        h1$k <- h1$Effect_D/h1$Effect_A
        h1$snpid <- gsub("S", "", h1$snpid)
        
        myt <- gsub(".*\\/|_.*", "", files[i])
        subtav <- subset(tav, trait == myt)
        
        subh <- merge(subtav, h1, by="snpid")
        message(sprintf("###>>> processing file: [ %s ] and get [ %s ] SNPs", files[i], nrow(subh)))
        
        tem <- subh$k
        out <- data.frame(file=files[i], 
                          k1=length(tem[tem > 1]),
                          k2=length(tem[tem >= 0.5 & tem <= 1]), 
                          k3=length(tem[tem > -0.5 & tem < 0.5]),
                          k4=length(tem[tem >= -1 & tem <= -0.5]), 
                          k5=length(tem[tem < -1]),
                          h2_A=sum(subh$h2_mrk_A),
                          h2_D=sum(subh$h2_mrk_D),
                          h2_tot=sum(subh$H2_mrk))
        
        output1 <- rbind(output1, out)
        output2 <- rbind(output2, subh)
        
    }
    return(list(output1, output2))
}


files <- list.files(path="largedata/snpeff", pattern="snpe$", full.names=TRUE)
tav <- read.csv("~/Documents/Heterosis_GWAS/HGWAS_proj/reports/StableN.tav758.csv")
res <- get_AD_eff(files, tav)

write.table(res[[1]], "cache/TAV_add_dom_num.csv", sep=",", row.names=FALSE, quote=FALSE)
write.table(res[[2]], "cache/TAV_add_dom_eff.csv", sep=",", row.names=FALSE, quote=FALSE)


