### Jinliang Yang
### 

tav <- read.csv("cache/TAV_filtered_4plot.csv")

tav$k <- tav$Effect_D/tav$Effect_A2

t <- as.character(unique(tav$trait))

out <- data.frame()
for(ti in t){
    sub <- subset(tav, trait == ti)
    
    pos <- nrow(subset(sub, k > 0.5))
    neg <- nrow(subset(sub, k < -0.5))
    tem <- data.frame(trait=ti, pos=pos, neg=neg)
    out <- rbind(out, tem)
}


