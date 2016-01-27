### Jinliang Yang
### Nov 25th, 2014
### purpose: cv with recombination rate suggested by reviewer two

get_training_map <- function(){
    nammap <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/NAM_map_and_genos/NAM_map_20080419.csv")
    nammap <- nammap[, 1:5]
    phy <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/data/NAM_map_and_genos/NAM_1144SNPs_AGPv2_positions.csv")
    phy <- phy[, 1:6]
    
    map <- merge(nammap, phy, by.x="marker", by.y="SNP_NAME")
    map <- map[, c("marker", "ch", "cumulative", "AGPv2_pos")]
    names(map) <- c("Marker", "Chr", "Genetic", "Physical")
    map <- map[map$Physical != "unknown", ]
    map$Physical <- as.numeric(as.character(map$Physical))
    return(map)
}

map <- get_training_map()


##############################################
cvd_p2g <- function(map=map){
    source("lib/p2g.R")
    
    ### format data for predicting genetic position
    cvsnp <- read.csv("reports/Stable13_cv_summary.csv")
    cvsnp <- cvsnp[, c(1,3,4,2,5:18)]
    
    #Input file:  Three columns: 1.marker_name 2.chromosome (Integer) 
    #3. Physical position (Colname must be Physical);
    cvp <- p2g(predict=cvsnp, train=map) 
    cvsnp <- merge(cvp, cvsnp[, -2:-3], by.x="marker", by.y="snpid")
    return(cvsnp)
}

cvsnp <- cvd_p2g(map=map)

save(file="cache/cvsnp_p2g.RData", list=c("map", "cvsnp"))

