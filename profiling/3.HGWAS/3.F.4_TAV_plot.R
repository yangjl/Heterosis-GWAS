### Jinliang Yang
### Feb 27th, 2016

get_gene_action <- function(TAV, kcutoff=0.1){
    traits <- c("KRN", "CD", "AKW", "CL", "CW", "KC", "TKW")
    
    output <- data.frame()
    for(i in 1:7){
        out <- subset(TAV, trait %in% traits[i])
        tem <- out$k
        out <- data.frame(trait=traits[i], 
                          k1=length(tem[tem > 1]),
                          k2=length(tem[tem >= kcutoff & tem <= 1]), 
                          k3=length(tem[tem > -kcutoff & tem < kcutoff]),
                          k4=length(tem[tem >= -1 & tem <= -kcutoff]), 
                          k5=length(tem[tem < -1]),
                          h2_A=sum(out$h2_mrk_A),
                          h2_D=sum(out$h2_mrk_D),
                          h2_tot=sum(out$H2_mrk))
        output <- rbind(output, out)
    }
    return(output)
}

res <- get_gene_action(TAV, kcutoff=0.1)
    
output1 <- rbind(output1, out)
output2 <- rbind(output2, subh)


tav <- read.csv("~/Documents/Heterosis_GWAS/HGWAS_proj/reports/StableN.tav758.csv")
res <- get_AD_eff(files, tav)

write.table(res[[1]], "cache/TAV_add_dom_num.csv", sep=",", row.names=FALSE, quote=FALSE)
write.table(res[[2]], "cache/TAV_add_dom_eff.csv", sep=",", row.names=FALSE, quote=FALSE)
