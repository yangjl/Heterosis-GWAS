### Jinliang Yang
### Feb 27th, 2016
### CV results

library(ggplot2)

st1 <- read.csv("cache/diallel_cv_strategy1.csv")
st1$trait <- toupper(st1$trait)
st2 <- read.csv("cache/diallel_cv_strategy2.csv")
st2$trait <- toupper(st2$trait)

theme_set(theme_grey(base_size = 18)) 
p1 <- ggplot(st1, aes(factor(trait, c("KRN", "CD", "AKW", "CL", "CW", "KC", "TKW")), r)) + 
    geom_boxplot(aes(fill = factor(mode, labels=c("Additive", "Dominance")))) +
    ylim(0, 1) +
    xlab("") +
    ylab("Prediction Accuracy") +
    ggtitle("") + theme_bw() +
    labs(fill="Model") +
    theme(axis.text.x = element_text(angle = 0, size=15),
          axis.title=element_text(size=15,face="bold"))

p2 <- ggplot(st2, aes(factor(trait, c("KRN", "CD", "AKW", "CL", "CW", "KC", "TKW")), r)) + 
    geom_boxplot(aes(fill = factor(mode, labels=c("Additive", "Dominance")))) +
    ylim(0, 1) +
    xlab("") +
    ylab("Prediction Accuracy") +
    ggtitle("") + theme_bw() +
    labs(fill="Model") +
    theme(axis.text.x = element_text(angle = 0, size=15),
          axis.title=element_text(size=15,face="bold"))

########
source("~/Documents/Github/zmSNPtools/Rcodes/multiplot.R")
pdf("graphs/Fig3_CV_boxplot.pdf", width=12, height=5)
multiplot(p1, p2, cols=2)
dev.off()



###################################################################################

runttest <- function(res){
    
    myt <- c("KRN", "CD", "AKW", "CL", "CW", "KC", "TKW")
    output <- data.frame()
    for(i in 1:7){
        add <- subset(res, mode == "add" & trait == myt[i])
        dom <- subset(res, mode == "dom" & trait == myt[i])
        
        test <- t.test(dom$r, add$r)
        
        out <- data.frame(trait=myt[i], mdom=mean(dom$r), madd=mean(add$r), pval=test$p.value)
        output <- rbind(output, out)
    }
    output$diff <- output$mdom - output$madd
    output$qval <- p.adjust(output$pval, method="bonferroni")
    return(output)
}


res1 <- runttest(st1)
write.table(res1, "cache/cv_st1_pvals.csv", sep=",", row.names=FALSE, quote=FALSE)

res2 <- runttest(st2)
write.table(res2, "cache/cv_st2_pvals.csv", sep=",", row.names=FALSE, quote=FALSE)

