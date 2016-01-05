### Jinliang Yang
### select the trait-associated variants using a thinning procedure
### 8/15/2014

ob1 <- load("~/Documents/Heterosis_GWAS/HGWAS_proj/cache/snptestpval.RData")
ob2 <- load("~/Documents/Heterosis_GWAS/HGWAS_proj/cache/bayes.RData")
ob3 <- load("~/Documents/Heterosis_GWAS/HGWAS_proj/cache/swpval.RData")

source("~/Documents/Heterosis_GWAS/HGWAS_proj/lib/TAVselect.R")



cutoff1=0.1
cutoff2=50
cutoff3=10
topnum =3

lstkw <- TAVselect(snptest=tkw, bayes=tkw2, sw=tkw3, 
                   bayescutoff=cutoff1, snptestcutoff=cutoff2, cwcutoff=cutoff3, 
                   binsize=1e7, topinbin=topnum)
# input: SNPTEST [12966279], Bayes [1296628] and SW [300]
# >cutoff: SNPTEST [21], Bayes [7] and SW [16]
# after thinning: SNPTEST [21], Bayes [7] and SW [6]

lsakw <- TAVselect(snptest=akw, bayes=akw2, sw=akw3, 
                   bayescutoff=cutoff1, snptestcutoff=cutoff2, cwcutoff=cutoff3, 
                   binsize=1e7, topinbin=topnum)
# input: SNPTEST [12966279], Bayes [1296628] and SW [300]
# >cutoff: SNPTEST [18], Bayes [35] and SW [44]
# after thinning: SNPTEST [18], Bayes [35] and SW [43]
lskc <- TAVselect(snptest=kc, bayes=kc2, sw=kc3, 
                   bayescutoff=cutoff1, snptestcutoff=cutoff2, cwcutoff=cutoff3, 
                   binsize=1e7, topinbin=topnum)
# input: SNPTEST [12966279], Bayes [1296628] and SW [300]
# >cutoff: SNPTEST [52], Bayes [16] and SW [25]
# after thinning: SNPTEST [51], Bayes [16] and SW [25]
lscw <- TAVselect(snptest=cw, bayes=cw2, sw=cw3, 
                  bayescutoff=cutoff1, snptestcutoff=cutoff2, cwcutoff=cutoff3, 
                  binsize=1e7, topinbin=topnum)
# input: SNPTEST [12966279], Bayes [1296628] and SW [300]
# >cutoff: SNPTEST [5], Bayes [19] and SW [35]
# after thinning: SNPTEST [5], Bayes [19] and SW [26]
lscl <- TAVselect(snptest=cl, bayes=cl2, sw=cl3, 
                  bayescutoff=cutoff1, snptestcutoff=cutoff2, cwcutoff=cutoff3, 
                  binsize=1e7, topinbin=topnum)
# input: SNPTEST [12966279], Bayes [1296628] and SW [300]
# >cutoff: SNPTEST [55], Bayes [22] and SW [33]
# after thinning: SNPTEST [55], Bayes [22] and SW [33]
lscd <- TAVselect(snptest=cd, bayes=cd2, sw=cd3, 
                  bayescutoff=cutoff1, snptestcutoff=cutoff2, cwcutoff=cutoff3, 
                  binsize=1e7, topinbin=topnum)
# input: SNPTEST [12966279], Bayes [1296628] and SW [300]
# >cutoff: SNPTEST [2467], Bayes [55] and SW [43]
# after thinning: SNPTEST [78], Bayes [54] and SW [43]

krn <- read.csv("~/Documents/KRN_GWAS_v3/GWAS3_proj/reports/S.table5.KAV.bins_v2.csv")
names(krn)[1:3] <- c("snpid", "chr", "pos")
krn1 <- subset(krn, Method=="single-variant")
names(krn1)[4] <- "log10p"
krn2 <- subset(krn, Method=="bayes")
names(krn2)[4] <- "ModelFreq"
krn3 <- subset(krn, Method=="stepwise")
names(krn3)[4] <- "Pvalue"
lskrn <- TAVselect(snptest=krn1, bayes=krn2, sw=krn3, 
                   bayescutoff=cutoff1, snptestcutoff=cutoff2, cwcutoff=cutoff3, 
                   binsize=1e7, topinbin=topnum)
# input: SNPTEST [253], Bayes [433] and SW [300]
# >cutoff: SNPTEST [62], Bayes [142] and SW [63]
# after thinning: SNPTEST [37], Bayes [119] and SW [63]

source("~/Documents/Rcodes/save.append.R")
save.append(list=c("lskrn", "lscd", "lscw", "lscl", "lsakw", "lstkw", "lskc"), 
            file="~/Documents/Heterosis_GWAS/HGWAS_proj/cache/TAVs_7traits_v2.RData",
            description="use bayes=0.1, snptest=50, sw=10")

######
ob <- load(file="largedata/lcache/TAVs_7traits.RData")
head(lscd[[1]])






