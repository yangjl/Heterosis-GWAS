# Jinliang Yang
# Jan 12th, 2015
# harvest the results of model training with gerp and random SNPs


harvestCV <- function(dir="largedata/scripts/", fileptn="\\.ghatREL1"){
    
    files <- list.files(path = dir, pattern=fileptn, full.names=TRUE)
    ## file line of the shell file:
    message(sprintf("[ %s ] files detected!", length(files)))
    
    res <- unlist(lapply(1:length(files), function(i){
        #genotype gHat DTP Fix  meanBias PEV=Var(g/y)   R^2 
        ghat <- read.table(files[i], skip=1, header=FALSE)
        if(i %% 10000 ==0) {
            # Print on the screen some message
            message(sprintf("###>>> finished reading [ %s ] files!", i))
        }
        return(cor(ghat$V2, ghat$V3))
    }))
    
    resout <- data.frame(file=files, r=res)
    
    return(resout)
}

########
SplitName <- function(infile=resout){
    
    infile$file <- as.character(infile$file)
    infile$file <- gsub(".*/", "", infile$file)
    
    #infile$cs <- gsub("_.*", "", infile$file)
    infile$trait <- gsub("_.*", "", infile$file)
    
    
    infile$mode <- gsub("_ran.*", "", infile$file) 
    infile$mode <- gsub(".*_", "", infile$mode)
    
    infile$rand <- gsub("_cv.*", "", infile$file) 
    infile$rand <- gsub(".*_ran", "", infile$rand)
    
    infile$cv <- gsub("\\..*", "", infile$file) 
    infile$cv <- gsub(".*_cv", "", infile$cv)
    
    print(table(infile$trait))
    return(infile)
}

#### extract with real data
res1 <- harvestCV(dir="largedata/scripts/", fileptn="\\.ghatREL")
res1 <- SplitName(infile=res1) #885
print(table(res1$trait))
write.table(res1, "cache/diallel_cv_strategy1.csv", sep=",", row.names=FALSE, quote=FALSE)

res2 <- harvestCV(dir="largedata/scripts_cv2/", fileptn="\\.ghatREL")
res2 <- SplitName(infile=res2) #885
print(table(res2$trait))
write.table(res2, "cache/diallel_cv_strategy2.csv", sep=",", row.names=FALSE, quote=FALSE)





library(plyr, lib="~/bin/Rlib/")


res1 <- collect_res(dir="slurm-scripts/cv_b2/")
test <- ddply(res1, .(mode, trait, type), summarise,
              r = mean(r))


################################################################
#7*3*5*100*11 = [1] 115500

g2 <- collect_res(dir="slurm-scripts/cv_b2/")
write.table(g2, "cache/g2_k_perse.csv", sep=",", row.names=FALSE, quote=FALSE)

g0 <- collect_res(dir="slurm-scripts/cv_b0/")
test <- ddply(g0, .(mode, trait, type), summarise,
              r = mean(r))

write.table(g0, "cache/g0_k_perse.csv", sep=",", row.names=FALSE, quote=FALSE)


###############
bph0 <- collect_res(dir="largedata/SNP/bph_b0_cs/")
test <- ddply(bph0, .(mode, trait, type), summarise,
              r = mean(r))

write.table(bph0, "cache/g0_k_bph_2016.csv", sep=",", row.names=FALSE, quote=FALSE)

b0 <- read.csv("cache/g0_k_bph.csv")

test <- ddply(b0, .(mode, trait, type), summarise,
              r = mean(r))

res02 <- read.csv("cache/g0_k_bph.csv")
res2 <- ddply(bph0, .(type, trait, mode, sp), summarise,
              r = mean(r),
              m = median(r))
############
pdf("graphs/Fig3_BPH_3plots_test.pdf", height=4, width=12)
par(mfrow=c(1,3))

mybean(res2, mymode = "a2", ylim=c(0, 1), main="Additive", ylab="Cross-validation Accuracy")
mybean(res2, mymode = "d2", ylim=c(0, 1), main="Dominance", ylab="Cross-validation Accuracy")
mybean(res2, mymode = "h2", ylim=c(0, 1), main="Incomplete Dominance", ylab="Cross-validation Accuracy")

dev.off()