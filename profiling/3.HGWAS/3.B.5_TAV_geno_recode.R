### Jinliang Yang
### 8/15/2014

ob <- load("~/Documents/Heterosis_GWAS/HGWAS_proj/cache/TAVs_7traits.RData")
getsnpid <- function(ob=lstkw){
  
  names(ob[[3]])[1] <- "snpid"
  names(ob[[2]])[1] <- "snpid"
  snpid <- c(ob[[1]]$snpid, ob[[2]]$snpid, ob[[3]]$snpid)
  message(sprintf("Total: [%s]; Unique: [%s]", length(snpid), length(unique(snpid))))
  return(snpid)
}
############################
snpid1 <- getsnpid(ob=lstkw)
snpid2 <- getsnpid(ob=lsakw)
snpid3 <- getsnpid(ob=lskc)
snpid4 <- getsnpid(ob=lscd)
snpid5 <- getsnpid(ob=lscw)
snpid6 <- getsnpid(ob=lscl)
snpid7 <- getsnpid(ob=lskrn)

allsnpid <- c(snpid1, snpid2, snpid3, snpid4, snpid5, snpid6, snpid7);
allsnpid <- unique(allsnpid)
#[1] 5764
write.table(t(allsnpid), "~/Documents/Heterosis_GWAS/Method/SNPTEST/allsnpid.txt",
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

#####
qctool -g geno_13M_chr1.bgen -og chr1_subsetted.gen -incl-rsids allsnpid.txt
qctool -g geno_13M_chr2.bgen -og chr2_subsetted.gen -incl-rsids allsnpid.txt
qctool -g geno_13M_chr3.bgen -og chr3_subsetted.gen -incl-rsids allsnpid.txt
qctool -g geno_13M_chr4.bgen -og chr4_subsetted.gen -incl-rsids allsnpid.txt
qctool -g geno_13M_chr5.bgen -og chr5_subsetted.gen -incl-rsids allsnpid.txt
qctool -g geno_13M_chr6.bgen -og chr6_subsetted.gen -incl-rsids allsnpid.txt
qctool -g geno_13M_chr7.bgen -og chr7_subsetted.gen -incl-rsids allsnpid.txt
qctool -g geno_13M_chr8.bgen -og chr8_subsetted.gen -incl-rsids allsnpid.txt
qctool -g geno_13M_chr9.bgen -og chr9_subsetted.gen -incl-rsids allsnpid.txt
qctool -g geno_13M_chr10.bgen -og chr10_subsetted.gen -incl-rsids allsnpid.txt

#############################################
source("~/Documents/Heterosis_GWAS/HGWAS_proj/lib/genReCode.R")
gen10 <- read.table("~/Documents/Heterosis_GWAS/Method/SNPTEST/chr10_subsetted.gen", header=FALSE)
genReCode(gen=gen10, code="raw", 
          outfile="~/Documents/Heterosis_GWAS/Method/genoTAV/chr10_recoded.csv")

#############################################
source("~/Documents/Heterosis_GWAS/HGWAS_proj/lib/genReCode.R")
gen9 <- read.table("~/Documents/Heterosis_GWAS/Method/SNPTEST/chr9_subsetted.gen", header=FALSE)
genReCode(gen=gen9, code="raw", 
          outfile="~/Documents/Heterosis_GWAS/Method/genoTAV/chr9_recoded.csv")

gen8 <- read.table("~/Documents/Heterosis_GWAS/Method/SNPTEST/chr8_subsetted.gen", header=FALSE)
genReCode(gen=gen8, code="raw", 
          outfile="~/Documents/Heterosis_GWAS/Method/genoTAV/chr8_recoded.csv")

gen7 <- read.table("~/Documents/Heterosis_GWAS/Method/SNPTEST/chr7_subsetted.gen", header=FALSE)
genReCode(gen=gen7, code="raw", 
          outfile="~/Documents/Heterosis_GWAS/Method/genoTAV/chr7_recoded.csv")

gen6 <- read.table("~/Documents/Heterosis_GWAS/Method/SNPTEST/chr6_subsetted.gen", header=FALSE)
genReCode(gen=gen6, code="raw", 
          outfile="~/Documents/Heterosis_GWAS/Method/genoTAV/chr6_recoded.csv")
gen5 <- read.table("~/Documents/Heterosis_GWAS/Method/SNPTEST/chr5_subsetted.gen", header=FALSE)
genReCode(gen=gen5, code="raw", 
          outfile="~/Documents/Heterosis_GWAS/Method/genoTAV/chr5_recoded.csv")
#####
gen4 <- read.table("~/Documents/Heterosis_GWAS/Method/SNPTEST/chr4_subsetted.gen", header=FALSE)
genReCode(gen=gen4, code="raw", 
          outfile="~/Documents/Heterosis_GWAS/Method/genoTAV/chr4_recoded.csv")
gen3 <- read.table("~/Documents/Heterosis_GWAS/Method/SNPTEST/chr3_subsetted.gen", header=FALSE)
genReCode(gen=gen3, code="raw", 
          outfile="~/Documents/Heterosis_GWAS/Method/genoTAV/chr3_recoded.csv")
gen2 <- read.table("~/Documents/Heterosis_GWAS/Method/SNPTEST/chr2_subsetted.gen", header=FALSE)
genReCode(gen=gen2, code="raw", 
          outfile="~/Documents/Heterosis_GWAS/Method/genoTAV/chr2_recoded.csv")
gen1 <- read.table("~/Documents/Heterosis_GWAS/Method/SNPTEST/chr1_subsetted.gen", header=FALSE)
genReCode(gen=gen1, code="raw", 
          outfile="~/Documents/Heterosis_GWAS/Method/genoTAV/chr1_recoded.csv")












