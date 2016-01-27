### Jinliang Yang
### Nov 25th, 2014
### purpose: cv with recombination rate suggested by reviewer two

get_training_map <- function(){
    nammap <- read.csv("data/NAM_map_20080419.csv")
    nammap <- nammap[, 1:5]
    phy <- read.csv("data/NAM_1144SNPs_AGPv2_positions.csv")
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
    tav <- read.csv("reports/StableN.tav758.csv")
    
    #Input file:  Three columns: 1.marker_name 2.chromosome (Integer) 
    #3. Physical position (Colname must be Physical);
    cvp <- p2g(predict=tav, train=map) 
    cvsnp <- merge(cvp, tav[, -2:-3], by.x="marker", by.y="snpid")
    return(cvsnp)
}

cvsnp <- cvd_p2g(map=map)

save(file="cache/cvsnp_p2g_TAV758.RData", list=c("map", "cvsnp"))

