

#####
cd largedata/SNP/
qctool -g geno_13M_chr1.bgen -og chr1_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr2.bgen -og chr2_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr3.bgen -og chr3_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr4.bgen -og chr4_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr5.bgen -og chr5_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr6.bgen -og chr6_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr7.bgen -og chr7_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr8.bgen -og chr8_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr9.bgen -og chr9_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr10.bgen -og chr10_tavsnp.gen -incl-rsids TAV_snpid.txt

#############################################
source("lib/genReCode.R")

chr <- data.frame()
for(i in 1:10){
    gen <- read.table(paste0("largedata/SNP/chr", i, "_tavsnp.gen"), header=FALSE)
    message(sprintf("###>>> recode [ chr %s ]", i))
    out <- genReCode(gen=gen, code="code0123", outfile=NULL)
    chr <- rbind(chr, out)
}

write.table(chr, "largedata/SNP/TAV_recoded.csv", sep=",", row.names=FALSE, quote=FALSE)


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




