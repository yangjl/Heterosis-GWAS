### Jinliang Yang
### select the trait-associated variants using a thinning procedure
### 8/15/2014
### 1/4/2016

ob1 <- load("largedata/lcache/snptestpval.RData")
ob2 <- load("largedata/lcache/bayes.RData")
ob3 <- load("largedata/lcache/swpval.RData")

source("lib/TAVselect.R")

cutoff1=0.02
cutoff2=20
cutoff3=3
topnum =10


lstkw <- TAVselect(snptest=tkw, bayes=tkw2, sw=tkw3, 
                   bayescutoff=0.02, snptestcutoff=20, cwcutoff=3, 
                   binsize=1e7, topinbin=10)
# input: SNPTEST [12966279], Bayes [1296628] and SW [300]
# >cutoff: SNPTEST [72], Bayes [114] and SW [172]
# after thinning: SNPTEST [72], Bayes [114] and SW [161]
lsakw <- TAVselect(snptest=akw, bayes=akw2, sw=akw3, 
                   bayescutoff=0.02, snptestcutoff=20, cwcutoff=3, 
                   binsize=1e7, topinbin=10)
# input: SNPTEST [12966279], Bayes [1296628] and SW [300]
# >cutoff: SNPTEST [182], Bayes [410] and SW [241]
# after thinning: SNPTEST [182], Bayes [374] and SW [241]
lskc <- TAVselect(snptest=kc, bayes=kc2, sw=kc3, 
                   bayescutoff=0.02, snptestcutoff=20, cwcutoff=3, 
                   binsize=1e7, topinbin=10)
# input: SNPTEST [12966279], Bayes [1296628] and SW [300]
# >cutoff: SNPTEST [6766], Bayes [234] and SW [214]
# after thinning: SNPTEST [347], Bayes [225] and SW [214]
lscw <- TAVselect(snptest=cw, bayes=cw2, sw=cw3, 
                  bayescutoff=0.02, snptestcutoff=20, cwcutoff=3, 
                  binsize=1e7, topinbin=10)
# input: SNPTEST [12966279], Bayes [1296628] and SW [300]
# >cutoff: SNPTEST [15066], Bayes [346] and SW [223]
# after thinning: SNPTEST [353], Bayes [321] and SW [219]
lscl <- TAVselect(snptest=cl, bayes=cl2, sw=cl3, 
                  bayescutoff=0.02, snptestcutoff=20, cwcutoff=3, 
                  binsize=1e7, topinbin=10)
# input: SNPTEST [12966279], Bayes [1296628] and SW [300]
# >cutoff: SNPTEST [12972], Bayes [283] and SW [232]
# after thinning: SNPTEST [467], Bayes [274] and SW [232]
lscd <- TAVselect(snptest=cd, bayes=cd2, sw=cd3, 
                  bayescutoff=0.02, snptestcutoff=20, cwcutoff=3, 
                  binsize=1e7, topinbin=10)
# input: SNPTEST [12966279], Bayes [1296628] and SW [300]
# >cutoff: SNPTEST [72351], Bayes [656] and SW [262]
# after thinning: SNPTEST [564], Bayes [586] and SW [262]

krn <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table5.KAV.bins_v2.csv")
names(krn)[1:3] <- c("snpid", "chr", "pos")
krn1 <- subset(krn, Method=="single-variant")
names(krn1)[4] <- "log10p"
krn2 <- subset(krn, Method=="bayes")
names(krn2)[4] <- "ModelFreq"
krn3 <- subset(krn, Method=="stepwise")
names(krn3)[4] <- "Pvalue"
lskrn <- TAVselect(snptest=krn1, bayes=krn2, sw=krn3, 
                   bayescutoff=0.02, snptestcutoff=20, cwcutoff=3, 
                   binsize=1e7, topinbin=10)
# input: SNPTEST [253], Bayes [433] and SW [300]
# >cutoff: SNPTEST [164], Bayes [432] and SW [261]
# after thinning: SNPTEST [150], Bayes [427] and SW [261]

source("~/Documents/Rcodes/save.append.R")
save.append(list=c("lskrn", "lscd", "lscw", "lscl", "lsakw", "lstkw", "lskc"), 
            file="~/Documents/Heterosis_GWAS/HGWAS_proj/cache/TAVs_7traits.RData",
            description="")

######
ob <- load(file="largedata/lcache/TAVs_7traits.RData")
head(lscd[[1]])






